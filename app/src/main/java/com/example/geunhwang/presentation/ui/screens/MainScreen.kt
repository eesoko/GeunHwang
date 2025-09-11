package com.example.geunhwang.presentation.ui.screens

import androidx.compose.foundation.layout.*
import androidx.compose.runtime.Composable
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.res.painterResource
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import androidx.wear.compose.material.*
import com.example.geunhwang.R

@Composable
fun MainScreen(
    onStartWorkoutClick: () -> Unit,
    onLogbookClick: () -> Unit,
    onDataCollectionClick: () -> Unit // 개발용
) {
    Column(
        modifier = Modifier
            .fillMaxSize()
            .padding(horizontal = 16.dp),
        verticalArrangement = Arrangement.Center,
        horizontalAlignment = Alignment.CenterHorizontally
    ) {
        Text("근황", fontSize = 24.sp, fontWeight = FontWeight.Bold)
        Spacer(modifier = Modifier.height(16.dp))

        // 운동 시작 버튼
        Button(
            onClick = onStartWorkoutClick,
            modifier = Modifier.fillMaxWidth(),
            colors = ButtonDefaults.buttonColors(backgroundColor = MaterialTheme.colors.primary)
        ) {
            Row(verticalAlignment = Alignment.CenterVertically) {
                Icon(
                    painter = painterResource(id = R.drawable.ic_fitness),
                    contentDescription = "운동 시작"
                )
                Spacer(Modifier.width(8.dp))
                Text("운동 시작")
            }
        }

        Spacer(modifier = Modifier.height(8.dp))

        // 운동 기록 버튼
        Button(
            onClick = onLogbookClick,
            modifier = Modifier.fillMaxWidth(),
            colors = ButtonDefaults.buttonColors(backgroundColor = Color(0xFF0056B3))
        ) {
            Row(verticalAlignment = Alignment.CenterVertically) {
                Icon(
                    painter = painterResource(id = R.drawable.ic_history),
                    contentDescription = "운동 기록"
                )
                Spacer(Modifier.width(8.dp))
                Text("운동 기록")
            }
        }

        Spacer(modifier = Modifier.height(8.dp))

        // (개발용) 데이터 수집 버튼
        Button(
            onClick = onDataCollectionClick,
            modifier = Modifier.fillMaxWidth(),
            colors = ButtonDefaults.buttonColors(backgroundColor = Color.DarkGray)
        ) {
            Row(verticalAlignment = Alignment.CenterVertically) {
                Icon(
                    painter = painterResource(id = R.drawable.ic_data_collect),
                    contentDescription = "데이터 수집"
                )
                Spacer(Modifier.width(8.dp))
                Text("데이터 수집")
            }
        }
    }
}