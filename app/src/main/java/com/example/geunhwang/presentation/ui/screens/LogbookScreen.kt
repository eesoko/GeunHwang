package com.example.geunhwang.presentation.ui.screens

import androidx.compose.foundation.layout.*
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import androidx.wear.compose.foundation.lazy.ScalingLazyColumn
import androidx.wear.compose.foundation.lazy.items
import androidx.wear.compose.material.Text
import com.example.geunhwang.presentation.ui.LogbookEntry
import java.text.SimpleDateFormat
import java.util.*
import java.util.concurrent.TimeUnit

@Composable
fun LogbookScreen(
    // ViewModel로부터 가공된 운동 기록 리스트를 받습니다.
    entries: List<LogbookEntry>
) {
    Column(modifier = Modifier
        .fillMaxSize()
        .padding(16.dp)) {

        Text(
            text = SimpleDateFormat("yyyy.MM.dd EEE", Locale.getDefault()).format(Date()),
            fontWeight = FontWeight.Bold,
            modifier = Modifier
                .fillMaxWidth()
                .padding(bottom = 16.dp)
        )

        // 스크롤 가능한 운동 기록 목록
        ScalingLazyColumn(
            modifier = Modifier.fillMaxSize(),
            verticalArrangement = Arrangement.spacedBy(16.dp)
        ) {
            if (entries.isEmpty()) {
                item {
                    Text("오늘의 운동 기록이 없습니다.")
                }
            } else {
                // ViewModel에서 받아온 entries 리스트를 반복하며 각 항목을 표시합니다.
                items(entries) { entry ->
                    Column {
                        val totalTimeMinutes = TimeUnit.MILLISECONDS.toMinutes(entry.totalTimeMillis)
                        val totalTimeSeconds = TimeUnit.MILLISECONDS.toSeconds(entry.totalTimeMillis) % 60
                        val timeString = String.format("%02d:%02d", totalTimeMinutes, totalTimeSeconds)

                        // 헤더 라인
                        Text(
                            "${entry.exerciseName} ${entry.totalSets}Sets • $timeString",
                            fontWeight = FontWeight.Bold
                        )
                        Spacer(modifier = Modifier.height(4.dp))
                        // 데이터 라인
                        Text(
                            entry.setsDetail,
                            fontSize = 12.sp
                        )
                    }
                }
            }
        }
    }
}
