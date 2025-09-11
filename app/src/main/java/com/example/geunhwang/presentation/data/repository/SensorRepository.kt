package com.example.geunhwang.presentation.data.repository

import android.content.Context
import android.hardware.Sensor
import android.hardware.SensorEvent
import android.hardware.SensorEventListener
import android.hardware.SensorManager
import android.util.Log
import kotlinx.coroutines.*
import kotlin.math.sqrt

class SensorRepository(context: Context) : SensorEventListener {

    var onMotionStarted: (() -> Unit)? = null
    var onMotionStopped: ((List<FloatArray>) -> Unit)? = null

    private var sensorManager: SensorManager =
        context.getSystemService(Context.SENSOR_SERVICE) as SensorManager

    private var accelerometer: Sensor? =
        sensorManager.getDefaultSensor(Sensor.TYPE_LINEAR_ACCELERATION)

    // --- ▼▼▼ '규칙적인 펄스' 감지를 위한 새로운 알고리즘 변수들 ▼▼▼ ---
    private val PULSE_DETECTION_THRESHOLD = 16.0f // 1. 펄스 강도: 16.0f 이상의 움직임만 감지
    private val PULSE_COUNT_FOR_START = 2        // 2. 펄스 횟수: 2회 이상의 펄스가 감지되어야 운동 시작
    private val PULSE_WINDOW_DURATION_MS = 5000L // 3. 시간 창: 5초 이내에 발생해야 함
    private val SET_END_DURATION_MS = 3000L      // 4. 세트 종료: 3초 이상 펄스가 없으면 휴식으로 전환

    private var isListening = false
    private var isWorkoutMode = false // isMotionDetected -> isWorkoutMode 로 이름 변경
    private var motionStopJob: Job? = null
    private val detectionScope = CoroutineScope(Dispatchers.Default)

    private var sensorDataBuffer = mutableListOf<FloatArray>()
    private var pulseTimestamps = mutableListOf<Long>() // 펄스 발생 시간을 기록할 리스트
    private var isAboveThreshold = false // 펄스의 '정점'을 한 번만 감지하기 위한 플래그
    // --- ▲▲▲ ---

    fun startListening() {
        if (isListening) return
        isListening = true
        pulseTimestamps.clear()
        accelerometer?.let {
            sensorManager.registerListener(this, it, SensorManager.SENSOR_DELAY_GAME)
        }
    }

    fun stopListening() {
        if (!isListening) return
        isListening = false
        sensorManager.unregisterListener(this)
        motionStopJob?.cancel()
        sensorDataBuffer.clear()
    }

    override fun onSensorChanged(event: SensorEvent?) {
        if (event?.sensor?.type != Sensor.TYPE_LINEAR_ACCELERATION) return

        val values = event.values.clone()
        val magnitude = sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2])
        val currentTime = System.currentTimeMillis()

        Log.d("SensorRepo", "Magnitude: $magnitude")

        // 1. 움직임이 펄스 임계값을 넘었을 때
        if (magnitude > PULSE_DETECTION_THRESHOLD) {
            motionStopJob?.cancel() // '세트 종료' 타이머가 돌고 있었다면 취소 (아직 안 끝남)

            // 이전에 임계값 아래에 있다가 막 위로 올라온 '순간'을 포착 (펄스의 정점)
            if (!isAboveThreshold) {
                isAboveThreshold = true
                pulseTimestamps.add(currentTime)
                Log.d("SensorRepo", "Pulse detected! Count: ${pulseTimestamps.size}")
            }

            // 아직 운동 모드가 아니라면, 운동 시작 조건을 만족하는지 체크
            if (!isWorkoutMode) {
                // 오래된 타임스탬프 제거 (5초 이전 기록은 삭제)
                pulseTimestamps.removeAll { it < currentTime - PULSE_WINDOW_DURATION_MS }

                // 5초 안에 2개 이상의 펄스가 감지되었다면 -> 운동 시작!
                if (pulseTimestamps.size >= PULSE_COUNT_FOR_START) {
                    isWorkoutMode = true
                    sensorDataBuffer.clear()
                    onMotionStarted?.invoke()
                    Log.d("SensorRepo", "====== WORKOUT MODE STARTED ======")
                }
            }

            // 운동 모드일 때는 항상 데이터 버퍼에 기록
            if (isWorkoutMode) {
                sensorDataBuffer.add(values)
            }
        }
        // 2. 움직임이 펄스 임계값 아래로 내려왔을 때
        else {
            isAboveThreshold = false // 다시 펄스를 감지할 수 있도록 플래그 초기화

            // 운동 모드 중에 움직임이 멈췄다면 -> '세트 종료' 타이머 시작
            if (isWorkoutMode && motionStopJob?.isActive != true) {
                motionStopJob = detectionScope.launch {
                    delay(SET_END_DURATION_MS)
                    // 3초의 지연 후에도 계속 멈춰있다면 -> 진짜 세트 종료!
                    isWorkoutMode = false
                    pulseTimestamps.clear() // 다음 세트를 위해 펄스 카운트 초기화
                    onMotionStopped?.invoke(sensorDataBuffer.toList())
                    Log.d("SensorRepo", "====== WORKOUT MODE STOPPED (Set End) ======")
                }
            }
        }
    }

    override fun onAccuracyChanged(sensor: Sensor?, accuracy: Int) {}
}

