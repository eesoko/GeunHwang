package com.example.geunhwang.presentation.ui.screens

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Spacer
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.height
import androidx.compose.runtime.Composable
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import androidx.wear.compose.foundation.lazy.ScalingLazyColumn
import androidx.wear.compose.material.Chip
import androidx.wear.compose.material.ChipDefaults
import androidx.wear.compose.material.MaterialTheme
import androidx.wear.compose.material.Text

@Composable
fun LoggingScreen(exerciseName: String, onStopLoggingClick: () -> Unit) {
    ScalingLazyColumn(
        modifier = Modifier.fillMaxSize(),
        horizontalAlignment = Alignment.CenterHorizontally,
        verticalArrangement = Arrangement.Center
    ) {
        item {
            Text(text = exerciseName, fontSize = 24.sp)
        }
        item {
            Text(text = "데이터 기록 중...", fontSize = 20.sp, color = Color.Green)
        }
        item {
            Spacer(modifier = Modifier.height(20.dp))
        }
        item {
            // --- ▼▼▼ Button을 Chip으로 교체합니다 ▼▼▼ ---
            Chip(
                onClick = onStopLoggingClick,
                label = { Text("기록 중지") },
                modifier = Modifier.fillMaxWidth(0.8f), // 너비를 살짝 줄여 보기 좋게
                colors = ChipDefaults.chipColors(backgroundColor = MaterialTheme.colors.error)
            )
        }
    }
}

