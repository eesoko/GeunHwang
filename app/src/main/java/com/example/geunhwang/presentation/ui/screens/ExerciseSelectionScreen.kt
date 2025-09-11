package com.example.geunhwang.presentation.ui.screens

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.PaddingValues
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.runtime.Composable
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.wear.compose.foundation.lazy.ScalingLazyColumn
import androidx.wear.compose.foundation.lazy.rememberScalingLazyListState
import androidx.wear.compose.material.Chip
import androidx.wear.compose.material.ChipDefaults
import androidx.wear.compose.material.Text
import com.example.geunhwang.presentation.core.ExerciseType

@Composable
fun ExerciseSelectionScreen(onExerciseSelected: (String) -> Unit) {
    val listState = rememberScalingLazyListState()

    ScalingLazyColumn(
        state = listState,
        modifier = Modifier.fillMaxWidth(),
        horizontalAlignment = Alignment.CenterHorizontally,
        verticalArrangement = Arrangement.spacedBy(8.dp),
        contentPadding = PaddingValues(top = 40.dp, bottom = 40.dp, start = 8.dp, end = 8.dp)
    ) {
        item { Text(text = "수집할 운동 선택", modifier = Modifier.padding(bottom = 8.dp)) }

        // --- ▼▼▼ Button을 Chip으로 교체합니다 ▼▼▼ ---
        item { ExerciseChip(exerciseName = ExerciseType.SQUAT, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.OVERHEAD_PRESS, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.PUSH_UP, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.SIDE_LATERAL_RAISE, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.DUMBBELL_CURL, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.LUNGE, onClick = onExerciseSelected) }
        item { ExerciseChip(exerciseName = ExerciseType.DUMBBELL_ROW, onClick = onExerciseSelected) }
    }
}

// 재사용을 위해 별도의 Composable로 분리
@Composable
fun ExerciseChip(exerciseName: String, onClick: (String) -> Unit) {
    Chip(
        onClick = { onClick(exerciseName) },
        label = {
            Text(
                text = exerciseName,
                textAlign = TextAlign.Center, // 텍스트를 중앙 정렬
                modifier = Modifier.fillMaxWidth()
            )
        },
        modifier = Modifier.fillMaxWidth(),
        colors = ChipDefaults.primaryChipColors() // 기본 칩 스타일
    )
}

