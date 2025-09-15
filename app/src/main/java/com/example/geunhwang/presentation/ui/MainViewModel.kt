package com.example.geunhwang.presentation.ui

import android.app.Application
import android.content.Context
import android.hardware.Sensor
import android.hardware.SensorEvent
import android.hardware.SensorEventListener
import android.hardware.SensorManager
import androidx.lifecycle.AndroidViewModel
import androidx.lifecycle.viewModelScope
import com.example.geunhwang.presentation.data.SensorDataLogger
import com.example.geunhwang.presentation.data.local.AppDatabase
import com.example.geunhwang.presentation.data.local.WorkoutSet
import com.example.geunhwang.presentation.data.repository.SensorRepository
import com.example.geunhwang.presentation.data.repository.WorkoutRepository
import kotlinx.coroutines.Job
import kotlinx.coroutines.delay
import kotlinx.coroutines.flow.*
import kotlinx.coroutines.launch
import java.util.Calendar


// --- ▼▼▼ [새로운 클래스] C++ 분류기를 호출하는 통역사 클래스 ▼▼▼ ---
class ExerciseClassifier {
    // Companion object는 클래스의 모든 인스턴스가 공유하는 정적(static) 블록입니다.
    companion object {
        // 앱이 시작될 때 단 한 번만 C++ 라이브러리를 메모리에 로드합니다.
        // "geunhwang"은 CMakeLists.txt에서 add_library()로 정의한 이름과 일치해야 합니다.
        init {
            System.loadLibrary("geunhwang")
        }
    }

    /**
     * 이 함수는 C++에 구현된 네이티브 함수를 호출하는 '다리' 역할을 합니다.
     * 'external' 키워드는 "이 함수의 실제 내용은 Kotlin이 아닌 외부에 있다"고 알려줍니다.
     * 함수 시그니처(이름, 파라미터, 반환 타입)는 C++의 JNI 함수와 정확히 일치해야 합니다.
     * @param sensorData SensorRepository로부터 받은 2차원 센서 데이터 배열
     * @return C++ 함수가 반환하는 "운동이름,횟수" 형태의 문자열
     */
    external fun classify(sensorData: Array<FloatArray>): String
}


// ... (LogbookEntry, WorkoutState 정의는 기존과 동일) ...
data class LogbookEntry(
    val exerciseName: String,
    val totalSets: Int,
    val totalTimeMillis: Long,
    val setsDetail: String
)
sealed class WorkoutState {
    object InitialRest : WorkoutState()
    object WorkingOut : WorkoutState()
    data class MidWorkoutRest(
        val lastExercise: String,
        val lastWeight: Double,
        val repHistory: List<Int>
    ) : WorkoutState()
}

class MainViewModel(application: Application) : AndroidViewModel(application), SensorEventListener {

    // --- StateFlow 정의 (기존과 동일) ---
    private val _uiState = MutableStateFlow<WorkoutState>(WorkoutState.InitialRest)
    val uiState = _uiState.asStateFlow()
    private val _totalWorkoutTime = MutableStateFlow(0L)
    val totalWorkoutTime = _totalWorkoutTime.asStateFlow()
    private val _restTime = MutableStateFlow(0L)
    val restTime = _restTime.asStateFlow()
    private val _currentWeight = MutableStateFlow(20.0)
    val currentWeight = _currentWeight.asStateFlow()
    private val _logbookEntries = MutableStateFlow<List<LogbookEntry>>(emptyList())
    val logbookEntries = _logbookEntries.asStateFlow()

    // --- Repository 및 내부 변수 정의 (기존과 동일) ---
    private val workoutRepository: WorkoutRepository
    private val sensorRepository: SensorRepository
    private var totalTimerJob: Job? = null
    private var restTimerJob: Job? = null
    private var currentSessionSets = mutableMapOf<String, MutableList<Int>>()

    // --- ▼▼▼ AI 분류기 인스턴스 추가 ▼▼▼ ---
    private val exerciseClassifier: ExerciseClassifier

    private val sensorDataLogger: SensorDataLogger
    private val sensorManager: SensorManager
    private val accelerometer: Sensor?
    private val gyroscope: Sensor?
    private var lastAccelData = floatArrayOf(0f, 0f, 0f)
    private var lastGyroData = floatArrayOf(0f, 0f, 0f)

    init {
        val app = getApplication<Application>()
        val workoutDao = AppDatabase.getDatabase(app).workoutDao()
        workoutRepository = WorkoutRepository(workoutDao)
        sensorRepository = SensorRepository(app)

        // --- ▼▼▼ 분류기 인스턴스 생성 ▼▼▼ ---
        exerciseClassifier = ExerciseClassifier()

        sensorDataLogger = SensorDataLogger(app)
        sensorManager = app.getSystemService(Context.SENSOR_SERVICE) as SensorManager
        accelerometer = sensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER)
        gyroscope = sensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE)

        // --- ▼▼▼ [핵심 변경] onMotionStopped 콜백에서 AI 분류기 호출 ▼▼▼ ---
        sensorRepository.onMotionStopped = { sensorData ->
            // 1. C++ 함수를 호출하여 분석 결과를 받습니다. (List<FloatArray>를 Array<FloatArray>로 변환)
            val resultString = exerciseClassifier.classify(sensorData.toTypedArray())

            // 2. C++이 반환한 "운동이름,횟수" 문자열을 파싱합니다.
            val parts = resultString.split(",")
            val detectedExercise = parts.getOrNull(0)?.split(":")?.getOrNull(1)?.trim() ?: "Unknown"
            val detectedReps = parts.getOrNull(1)?.split(":")?.getOrNull(1)?.trim()?.toIntOrNull() ?: 0

            // 3. 분석된 결과로 세트 종료 로직을 처리합니다.
            onSetFinished(detectedExercise, detectedReps)
        }

        sensorRepository.onMotionStarted = { onWorkoutPulseDetected() }
        loadLogbookData()
    }

    // --- (이하 나머지 함수들은 기존과 동일) ---
    fun startDataLogging(exerciseName: String) {
        sensorDataLogger.startLogging(exerciseName)
        sensorManager.registerListener(this, accelerometer, SensorManager.SENSOR_DELAY_GAME)
        sensorManager.registerListener(this, gyroscope, SensorManager.SENSOR_DELAY_GAME)
    }

    fun stopDataLogging() {
        sensorManager.unregisterListener(this)
        sensorDataLogger.stopLogging()
    }

    override fun onAccuracyChanged(sensor: Sensor?, accuracy: Int) {}

    override fun onSensorChanged(event: SensorEvent?) {
        val timestamp = System.currentTimeMillis()
        when (event?.sensor?.type) {
            Sensor.TYPE_ACCELEROMETER -> lastAccelData = event.values.clone()
            Sensor.TYPE_GYROSCOPE -> lastGyroData = event.values.clone()
        }
        sensorDataLogger.logData(timestamp, lastAccelData, lastGyroData)
    }

    fun startWorkoutSession() { sensorRepository.startListening() }
    fun loadLogbookData() {
        viewModelScope.launch {
            workoutRepository.getTodaySets().collect { sets ->
                val groupedSets = sets.groupBy { it.exerciseName }

                val entries = groupedSets.map { (name, setsForExercise) ->
                    val sortedSets = setsForExercise.sortedBy { it.timestamp }
                    val totalTime = if (sortedSets.size > 1) {
                        sortedSets.last().timestamp - sortedSets.first().timestamp
                    } else 0

                    LogbookEntry(
                        exerciseName = name,
                        totalSets = setsForExercise.size,
                        totalTimeMillis = totalTime,
                        setsDetail = sortedSets.joinToString(" / ") { "${it.reps} ${it.weight}kg" }
                    )
                }
                _logbookEntries.value = entries.sortedByDescending { it.exerciseName }
            }
        }
    }
    fun onWorkoutPulseDetected() {
        if (_uiState.value is WorkoutState.WorkingOut) return
        _uiState.value = WorkoutState.WorkingOut
        stopRestTimer()
        startTotalTimer()
    }
    fun onSetFinished(exerciseName: String, reps: Int) {
        viewModelScope.launch {
            val weight = _currentWeight.value
            workoutRepository.insertSet(WorkoutSet(exerciseName = exerciseName, reps = reps, weight = weight))
            val history = currentSessionSets.getOrPut(exerciseName) { mutableListOf() }
            history.add(reps)
            _uiState.value = WorkoutState.MidWorkoutRest(
                lastExercise = exerciseName,
                lastWeight = weight,
                repHistory = history.toList()
            )
            startRestTimer()
        }
    }
    fun onWeightChange(delta: Double) { _currentWeight.update { (it + delta).coerceAtLeast(0.0) } }
    fun onWorkoutEnd() {
        sensorRepository.stopListening()
        _uiState.value = WorkoutState.InitialRest
        stopTotalTimer()
        stopRestTimer()
        currentSessionSets.clear()
        _totalWorkoutTime.value = 0L
        _restTime.value = 0L
        _currentWeight.value = 20.0
        loadLogbookData()
    }
    private fun startTotalTimer() {
        if (totalTimerJob?.isActive == true) return
        totalTimerJob = viewModelScope.launch { while (true) { delay(1000); _totalWorkoutTime.value++ } }
    }
    private fun stopTotalTimer() { totalTimerJob?.cancel(); totalTimerJob = null }
    private fun startRestTimer() {
        stopRestTimer()
        restTimerJob = viewModelScope.launch { _restTime.value = 0; while (true) { delay(1000); _restTime.value++ } }
    }
    private fun stopRestTimer() { restTimerJob?.cancel(); restTimerJob = null }
    override fun onCleared() {
        super.onCleared()
        sensorManager.unregisterListener(this)
        sensorRepository.stopListening()
        stopTotalTimer()
        stopRestTimer()
    }
}