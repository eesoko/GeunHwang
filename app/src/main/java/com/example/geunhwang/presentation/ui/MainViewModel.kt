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

// ... (LogbookEntry, WorkoutState 정의는 동일) ...
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

    // --- ▼▼▼ 데이터 수집을 위한 변수 추가 ▼▼▼ ---
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

        // --- 데이터 수집용 센서 및 로거 초기화 ---
        sensorDataLogger = SensorDataLogger(app)
        sensorManager = app.getSystemService(Context.SENSOR_SERVICE) as SensorManager
        accelerometer = sensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER)
        gyroscope = sensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE)

        // --- 센서-ViewModel 연결 (기존과 동일) ---
        sensorRepository.onMotionStarted = { onWorkoutPulseDetected() }
        sensorRepository.onMotionStopped = { sensorData ->
            val detectedExercise = "Squat"
            val detectedReps = 12
            onSetFinished(detectedExercise, detectedReps)
        }

        loadLogbookData()
    }

    // --- ▼▼▼ 데이터 수집을 위한 함수들 추가 ▼▼▼ ---
    fun startDataLogging(exerciseName: String) {
        sensorDataLogger.startLogging(exerciseName)
        sensorManager.registerListener(this, accelerometer, SensorManager.SENSOR_DELAY_GAME)
        sensorManager.registerListener(this, gyroscope, SensorManager.SENSOR_DELAY_GAME)
    }

    fun stopDataLogging() {
        sensorManager.unregisterListener(this)
        sensorDataLogger.stopLogging()
    }

    // --- SensorEventListener 인터페이스 구현 ---
    override fun onAccuracyChanged(sensor: Sensor?, accuracy: Int) {}

    override fun onSensorChanged(event: SensorEvent?) {
        // 이 함수는 데이터 수집 모드에서만 사용됩니다.
        val timestamp = System.currentTimeMillis()
        when (event?.sensor?.type) {
            Sensor.TYPE_ACCELEROMETER -> lastAccelData = event.values.clone()
            Sensor.TYPE_GYROSCOPE -> lastGyroData = event.values.clone()
        }
        // 양쪽 센서 값이 모두 준비되면 파일에 기록
        sensorDataLogger.logData(timestamp, lastAccelData, lastGyroData)
    }

    // --- (이하 나머지 함수들은 기존과 동일) ---
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