package com.example.geunhwang.presentation.ui

import android.os.Bundle
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.compose.runtime.*
import androidx.lifecycle.ViewModel
import androidx.lifecycle.ViewModelProvider
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
                // ViewModel을 생성할 때 MainActivity의 인스턴스(this)를 전달하는 Factory를 사용합니다.
                val viewModel: MainViewModel = viewModel(factory = MainViewModel.Factory(this))
                WearAppV2Navigation(viewModel = viewModel)
            }
        }
    }

    // --- 수정된 부분: internal -> private ---
    // C++ 네이티브 함수를 선언합니다. private으로 하여 이름 변경 문제를 방지합니다.
    // 이 함수는 ViewModel을 통해서만 호출될 것입니다.
    internal fun predictMotion(sensorData: Array<DoubleArray>, fs: Double): String {
        return predictMotionNative(sensorData, fs)
    }

    private external fun predictMotionNative(sensorData: Array<DoubleArray>, fs: Double): String

    // 네이티브 라이브러리를 로드합니다.
    companion object {
        init {
            System.loadLibrary("native-lib")
        }
    }
}


@Composable
fun WearAppV2Navigation(viewModel: MainViewModel) {
    val navController = rememberSwipeDismissableNavController()
    // ... (UI 상태 관찰 코드는 기존과 동일) ...
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
                onDataCollectionClick = { navController.navigate("collection_selection_screen") }
            )
        }

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

        composable("collection_selection_screen") {
            ExerciseSelectionScreen { exerciseName ->
                viewModel.startDataLogging(exerciseName)
                navController.navigate("logging_screen/$exerciseName")
            }
        }

        composable("logging_screen/{exerciseName}") { backStackEntry ->
            val exerciseName = backStackEntry.arguments?.getString("exerciseName") ?: "기록 중"
            LoggingScreen(
                exerciseName = exerciseName,
                onStopLoggingClick = {
                    viewModel.stopDataLogging()
                    navController.popBackStack()
                }
            )
        }
    }
}