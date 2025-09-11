package com.example.geunhwang.presentation.ui.screens

import androidx.compose.foundation.layout.*
import androidx.compose.runtime.Composable
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import androidx.wear.compose.material.Button
import androidx.wear.compose.material.ButtonDefaults
import androidx.wear.compose.material.Text
import androidx.wear.compose.material.TimeText
import com.example.geunhwang.presentation.ui.WorkoutState
import java.util.concurrent.TimeUnit

@Composable
fun WorkoutScreen(
    uiState: WorkoutState,
    totalWorkoutTime: Long,
    restTime: Long,
    currentWeight: Double,
    onWeightChange: (Double) -> Unit,
    onEndWorkoutClick: () -> Unit
) {
    Box(modifier = Modifier.fillMaxSize()) {
        TimeText(modifier = Modifier.align(Alignment.TopCenter))

        when (uiState) {
            is WorkoutState.InitialRest -> {
                InitialRestUI()
            }
            is WorkoutState.WorkingOut -> {
                WorkingOutUI()
            }
            is WorkoutState.MidWorkoutRest -> {
                MidWorkoutRestUI(
                    state = uiState,
                    totalWorkoutTime = totalWorkoutTime,
                    restTime = restTime,
                    currentWeight = currentWeight,
                    onWeightChange = onWeightChange
                )
            }
        }

        Button(
            onClick = onEndWorkoutClick,
            colors = ButtonDefaults.buttonColors(backgroundColor = Color.Red.copy(alpha = 0.8f)),
            modifier = Modifier
                .align(Alignment.BottomCenter)
                .fillMaxWidth()
                .padding(horizontal = 24.dp, vertical = 12.dp)
        ) {
            Text("운동 종료")
        }
    }
}

@Composable
fun InitialRestUI() {
    Box(modifier = Modifier.fillMaxSize(), contentAlignment = Alignment.Center) {
        Text("운동합시다!", fontSize = 24.sp, fontWeight = FontWeight.Bold)
    }
}

@Composable
fun WorkingOutUI() {
    Box(modifier = Modifier.fillMaxSize(), contentAlignment = Alignment.Center) {
        Text("분석중...", fontSize = 24.sp, fontWeight = FontWeight.Bold, color = Color.Green)
    }
}

@Composable
fun MidWorkoutRestUI(
    state: WorkoutState.MidWorkoutRest,
    totalWorkoutTime: Long,
    restTime: Long,
    currentWeight: Double,
    onWeightChange: (Double) -> Unit
) {
    Column(
        modifier = Modifier.fillMaxSize(),
        horizontalAlignment = Alignment.CenterHorizontally,
        verticalArrangement = Arrangement.Center
    ) {
        Row {
            Text("${restTime}s", color = Color.Yellow)
            Text(" • ")
            Text(formatTime(totalWorkoutTime))
        }
        Spacer(modifier = Modifier.height(12.dp))

        val lastReps = state.repHistory.lastOrNull() ?: 0
        Text("${state.lastExercise} ($lastReps reps)", fontSize = 22.sp, fontWeight = FontWeight.Bold)
        Spacer(modifier = Modifier.height(8.dp))

        Row(
            verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(8.dp)
        ) {
            Button(onClick = { onWeightChange(-1.0) }) { Text("-1") }
            Text("%.1f kg".format(currentWeight), fontSize = 20.sp)
            Button(onClick = { onWeightChange(1.0) }) { Text("+1") }
        }

        Spacer(modifier = Modifier.height(60.dp))
    }
}

private fun formatTime(seconds: Long): String {
    val minutes = TimeUnit.SECONDS.toMinutes(seconds)
    val remainingSeconds = seconds - TimeUnit.MINUTES.toSeconds(minutes)
    return String.format("%02d:%02d", minutes, remainingSeconds)
}

