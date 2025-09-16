package com.example.geunhwang.presentation.ui

import android.app.Application
import android.content.Context
import android.hardware.Sensor
import android.hardware.SensorEvent
import android.hardware.SensorEventListener
import android.hardware.SensorManager
import androidx.lifecycle.ViewModel
import androidx.lifecycle.ViewModelProvider
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

// ... (LogbookEntry, WorkoutState data/sealed class 정의는 기존과 동일) ...
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


// --- ▼▼▼ MainActivity 참조를 받도록 생성자 수정 ▼▼▼ ---
class MainViewModel(
    application: Application,
    private val mainActivity: MainActivity
) : ViewModel(), SensorEventListener {

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

    // --- 데이터 로깅 관련 변수 (기존과 동일) ---
    private val sensorDataLogger: SensorDataLogger
    private val sensorManager: SensorManager
    private val accelerometer: Sensor?
    private val gyroscope: Sensor?
    private var lastAccelData = floatArrayOf(0f, 0f, 0f)
    private var lastGyroData = floatArrayOf(0f, 0f, 0f)

    init {
        val app = application
        val workoutDao = AppDatabase.getDatabase(app).workoutDao()
        workoutRepository = WorkoutRepository(workoutDao)
        sensorRepository = SensorRepository(app)

        sensorDataLogger = SensorDataLogger(app)
        sensorManager = app.getSystemService(Context.SENSOR_SERVICE) as SensorManager
        accelerometer = sensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER)
        gyroscope = sensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE)

        // --- onMotionStopped 콜백에서 AI 분류기 호출 (기존과 동일) ---
        sensorRepository.onMotionStopped = { sensorData ->
            // 데이터가 비어있지 않을 때만 아래 로직을 실행하도록 변경
            if (sensorData.isNotEmpty()) {
                val doubleArray = sensorData.map { it.map { value -> value.toDouble() }.toDoubleArray() }.toTypedArray()
                val resultString = mainActivity.predictMotion(doubleArray, 50.0)

                val parts = resultString.split(",")
                val detectedExercise = parts.getOrNull(0) ?: "Unknown"
                val detectedReps = parts.getOrNull(1)?.toIntOrNull() ?: 0

                onSetFinished(detectedExercise, detectedReps)
            }
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

    // --- ▼▼▼ ViewModel에 Activity를 전달하기 위한 Factory 클래스 추가 ▼▼▼ ---
    class Factory(private val activity: MainActivity) : ViewModelProvider.Factory {
        override fun <T : ViewModel> create(modelClass: Class<T>): T {
            if (modelClass.isAssignableFrom(MainViewModel::class.java)) {
                @Suppress("UNCHECKED_CAST")
                return MainViewModel(activity.application, activity) as T
            }
            throw IllegalArgumentException("Unknown ViewModel class")
        }
    }
}