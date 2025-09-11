package com.example.geunhwang.presentation.ui

import android.os.Bundle
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.compose.runtime.*
import androidx.lifecycle.viewmodel.compose.viewModel
import androidx.wear.compose.navigation.SwipeDismissableNavHost
import androidx.wear.compose.navigation.composable
import androidx.wear.compose.navigation.rememberSwipeDismissableNavController
import com.example.geunhwang.presentation.ui.screens.ExerciseSelectionScreen
import com.example.geunhwang.presentation.ui.screens.LogbookScreen
import com.example.geunhwang.presentation.ui.screens.LoggingScreen
import com.example.geunhwang.presentation.ui.screens.MainScreen
import com.example.geunhwang.presentation.ui.screens.WorkoutScreen
import com.example.geunhwang.presentation.ui.theme.GeunHwangTheme

class MainActivity : ComponentActivity() {
    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContent {
            GeunHwangTheme {
                val viewModel: MainViewModel = viewModel()
                WearAppV2Navigation(viewModel = viewModel)
            }
        }
    }
}

@Composable
fun WearAppV2Navigation(viewModel: MainViewModel) {
    val navController = rememberSwipeDismissableNavController()
    // ... (state 관찰은 기존과 동일) ...
    val uiState by viewModel.uiState.collectAsState()
    val totalWorkoutTime by viewModel.totalWorkoutTime.collectAsState()
    val restTime by viewModel.restTime.collectAsState()
    val currentWeight by viewModel.currentWeight.collectAsState()
    val logbookEntries by viewModel.logbookEntries.collectAsState()

    SwipeDismissableNavHost(
        navController = navController,
        startDestination = "main_screen"
    ) {
        composable("main_screen") {
            MainScreen(
                onStartWorkoutClick = { navController.navigate("workout_screen") },
                onLogbookClick = {
                    viewModel.loadLogbookData()
                    navController.navigate("logbook_screen")
                },
                // --- ▼▼▼ 데이터 수집 버튼에 기능 연결 ▼▼▼ ---
                onDataCollectionClick = { navController.navigate("collection_selection_screen") }
            )
        }

        // --- (workout_screen, logbook_screen은 기존과 동일) ---
        composable("workout_screen") {
            DisposableEffect(Unit) {
                viewModel.startWorkoutSession()
                onDispose {
                    if (viewModel.uiState.value !is WorkoutState.InitialRest) {
                        viewModel.onWorkoutEnd()
                    }
                }
            }
            WorkoutScreen(
                uiState = uiState,
                totalWorkoutTime = totalWorkoutTime,
                restTime = restTime,
                currentWeight = currentWeight,
                onWeightChange = { delta -> viewModel.onWeightChange(delta) },
                onEndWorkoutClick = {
                    viewModel.onWorkoutEnd()
                    navController.popBackStack()
                }
            )
        }
        composable("logbook_screen") {
            LogbookScreen(entries = logbookEntries)
        }

        // --- ▼▼▼ 데이터 수집을 위한 새로운 화면 경로 추가 ▼▼▼ ---
        composable("collection_selection_screen") {
            ExerciseSelectionScreen { exerciseName ->
                // 운동 선택 시, ViewModel에 로깅 시작을 알리고 로깅 화면으로 이동
                viewModel.startDataLogging(exerciseName)
                navController.navigate("logging_screen/$exerciseName")
            }
        }

        composable("logging_screen/{exerciseName}") { backStackEntry ->
            val exerciseName = backStackEntry.arguments?.getString("exerciseName") ?: "기록 중"
            LoggingScreen(
                exerciseName = exerciseName,
                onStopLoggingClick = {
                    // 기록 중지 시, ViewModel에 로깅 중지를 알리고 이전 화면으로 복귀
                    viewModel.stopDataLogging()
                    navController.popBackStack()
                }
            )
        }
    }
}