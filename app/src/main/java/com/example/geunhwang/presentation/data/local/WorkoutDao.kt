package com.example.geunhwang.presentation.data.local

import androidx.room.Dao
import androidx.room.Insert
import androidx.room.Query
import kotlinx.coroutines.flow.Flow

@Dao
interface WorkoutDao {

    @Insert
    suspend fun insertSet(workoutSet: WorkoutSet)

    // 함수 이름을 원래의 getTodaySets로 되돌리고, ORDER BY DESC (최신순 정렬)만 추가합니다.
    @Query("SELECT * FROM workout_sets WHERE timestamp >= :todayTimestamp ORDER BY timestamp DESC")
    fun getTodaySets(todayTimestamp: Long): Flow<List<WorkoutSet>>
}

