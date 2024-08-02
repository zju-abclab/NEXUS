#include "Choosemax.cuh"

/**************************************************
 * find the index of minimum among c[0] ~ c[num-1]
 *************************************************/
int MinIndex(RR* c, int num) {
  int min_index = 0;
  RR min = c[0];
  for (int i = 1; i < num; i++) {
    if (c[i] < min) {
      min_index = i;
      min = c[i];
    }
  }
  return min_index;
}

/******************************************************
 * find n indices of alternating sign whose sum is maximum among a[0] ~ a[m-1]
 * after executing this function, cur_index[0] ~ cur_index[n-1] become the indices of maximum sum
 * ***************************************************/
void MaxSubsetSum(RR* a, int m, int n, int* cur_index) {
  int cur_m, min_index;
  RR* adj_sum = new RR[m];

  cur_m = m;
  for (int i = 0; i <= m - 1; i++) cur_index[i] = i;
  while (cur_m > n) {
    if (cur_m - n >= 3) {
      for (int i = 0; i <= cur_m - 2; i++) adj_sum[i] = a[cur_index[i]] + a[cur_index[i + 1]];
      min_index = MinIndex(adj_sum, cur_m - 1);
      if (min_index == 0) {
        for (int i = 0; i <= cur_m - 2; i++) cur_index[i] = cur_index[i + 1];
        cur_m -= 1;
      } else if (min_index == cur_m - 2) {
        cur_m -= 1;
      } else {
        for (int i = min_index; i <= cur_m - 3; i++)
          cur_index[i] = cur_index[i + 2];
        cur_m -= 2;
      }
    } else if (cur_m - n == 2) {
      for (int i = 0; i <= cur_m - 2; i++) adj_sum[i] = a[cur_index[i]] + a[cur_index[i + 1]];
      adj_sum[cur_m - 1] = a[cur_index[cur_m - 1]] + a[cur_index[0]];
      min_index = MinIndex(adj_sum, cur_m);
      if (min_index == cur_m - 1) {
        for (int i = 0; i <= cur_m - 3; i++) cur_index[i] = cur_index[i + 1];
        cur_m -= 2;
      } else {
        for (int i = min_index; i <= cur_m - 3; i++) cur_index[i] = cur_index[i + 2];
        cur_m -= 2;
      }
    } else if (cur_m - n == 1) {
      min_index = MinIndex(a, cur_m);
      for (int i = min_index; i <= cur_m - 2; i++) cur_index[i] = cur_index[i + 1];

      cur_m -= 1;
    }
  }
  delete[] adj_sum;
}
