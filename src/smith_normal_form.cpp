#include "SmithNormalForm.h"

int_t rank_, consistent, numEDs, choiceA, choiceB;
int_t maxED = 100;

void SNF(std::vector<std::vector<int_t>> &A, std::vector<std::vector<int_t>> &P,
         std::vector<std::vector<int_t>> &Q) {
  int_t N = A.size();
  int_t M = N;

  std::vector<std::vector<int_t>> Pinv(N, std::vector<int_t>(N, 0));
  std::vector<std::vector<int_t>> Qinv(N, std::vector<int_t>(N, 0));

  for (int_t i = 0; i < A.size(); ++i) {
    P[i][i] = 1;
    Pinv[i][i] = 1;
    Q[i][i] = 1;
    Qinv[i][i] = 1;
  }

  int_t n, m;                   // used in for loops
  int_t tempN, tempM, tempMin;  // used to store temporary locations in A
  int_t finished = 0;           // bool to decide if loop is over
  int_t divides, allEntriesPos, smith;  // 1 is true, 0 is false
  int_t tempRowEntry, tempColEntry, tempMult, tempEntry;
  int_t diag = 0;

  // for finishedRows and finishedColumns,
  // set entry i to -1 if row/col i is not finished, or
  // set entry i to i if row/col i is finished
  std::vector<int_t> finishedRows(N);
  std::vector<int_t> finishedColumns(M);

  // initializes finishedRows and finishedColumns
  for (n = 0; n < N; ++n) {
    finishedRows[n] = -1;
  }
  for (m = 0; m < M; ++m) {
    finishedColumns[m] = -1;
  }

  while (finished == 0) {
    // finds least non-zero entry, stores it in tempN, tempM
    tempM = -1;
    tempN = -1;
    findLeastEntry(A, finishedRows, finishedColumns, &tempN, &tempM, &finished);

    if (finished == 1) {
      break;
    }

    // checks if entry divides its entire row and column
    divides = dividesRowAndCol(A, tempN, tempM);

    /* if entry divides its entire row and column,
    we annihilate the entire row and column with euclidean division */
    if (divides == 1) {
      for (n = 0; n < N; ++n) {
        // we do not operate on the row itself
        if (n != tempN) {
          tempMult = eucDiv(A[n][tempM], A[tempN][tempM]);

          /* if tempMult is 0, this operation is useless
             so we do not execute it.*/
          if (tempMult != 0) {
            rowOperations2(A, P, Pinv, n, tempN, tempMult);
          }
        }
      }
      for (m = 0; m < M; ++m) {
        // we do not operate on the column itself
        if (m != tempM) {
          tempMult = eucDiv(A[tempN][m], A[tempN][tempM]);

          /* if tempMult is 0, this operation is useless
          so we do not execute it. */
          if (tempMult != 0) {
            columnOperations2(A, Q, Qinv, m, tempM, tempMult);
          }
        }
      }
    } else {
      tempRowEntry = -1;
      tempColEntry = -1;
      if (contains(finishedColumns, M, tempM) == 0) {
        for (n = 0; n < N; ++n) {
          if (A[n][tempM] != 0 && n != tempN) {
            tempRowEntry = n;
            break;
          }
        }
        if (tempRowEntry == -1) {
          break;
        }
        tempMult = eucDiv(A[tempRowEntry][tempM], A[tempN][tempM]);
        rowOperations2(A, P, Pinv, tempRowEntry, tempN, tempMult);
      } else {
        for (m = 0; m < M; ++m) {
          if (A[tempN][m] != 0 && m != tempM) {
            tempColEntry = m;
            break;
          }
        }
        if (tempColEntry == -1) {
          break;
        }
        tempMult = eucDiv(A[tempN][tempColEntry], A[tempN][tempM]);
        columnOperations2(A, Q, Qinv, tempColEntry, tempM, tempMult);
      }
    }

    // updates rows and colums that are finished
    updateFinishedRows(A, finishedRows);
    updateFinishedColumns(A, finishedColumns);
    finished = done(finishedRows, finishedColumns, N, M);
  }

  // rank_ = getRank(A);
  rank_ = A.size();

  /* if any entries are negative at this point_t,
     multiply whole row by -1 */
  makeAllDiagsPositive(A, P, Pinv);

  /* perform type 3 operations to order the
     non-zero elements onto the diagonals
  */
  orderDiagonals(A, P, Pinv, Q, Qinv);

  smith = checkSmith(A);
  while (smith != -1) {
    /* write something that takes a_n, a_n+1 and transforms
       it int_to gcd(a_n, a_n+1) and lcm(a_n, a_n+1) */
    smithTransform(A, P, Pinv, Q, Qinv, smith);
    // this call ensures all entries are ordered
    orderDiagonals(A, P, Pinv, Q, Qinv);
    smith = checkSmith(A);
  }
  // once again, ensure all diagonal elements are positive
  makeAllDiagsPositive(A, P, Pinv);

  // computation trick: transpose Q and Pinv so they are correct
  transposeM(Qinv);
  transposeN(Pinv);
}

void initializeZero(std::vector<int_t> &arr, int_t len) {
  int_t i;
  for (i = 0; i < len; ++i) {
    arr[i] = 0;
  }
}

// puts all the diagonal elements int_to an array
void getDiags(std::vector<std::vector<int_t>> &A, std::vector<int_t> &diags) {
  int_t i;
  for (i = 0; i < rank_; ++i) {
    diags[i] = A[i][i];
  }
}

// gets the elementary divisors using Prof. Zhou's algo
void getElementaryDivisors(std::vector<int_t> &diags, std::vector<int_t> &EDs) {
  int_t i, j;
  int_t ret;
  for (i = 0; i < rank_; ++i) {
    primeFactorize(diags[i], EDs);
  }
}

// returns absolute value of an int_teger
int_t absolu(int_t x) {
  if (x < 0) {
    x = -1 * x;
  }
  return x;
}

int_t min(int_t a, int_t b) {
  if (a < b) {
    return a;
  } else {
    return b;
  }
}

/*
to find inverses, we will use the rule that transpose(transpose(invP)) = invP
*/

// row operations

// multiplies entire row by unit for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type1rowM(std::vector<std::vector<int_t>> &A, int_t row, int_t unit) {
  int_t M = A.size();
  if (unit != 1 && unit != -1) {
    printf(
        "A type 1 row operation did not execute because the input was not a "
        "unit.");
    return;
  } else {
    int_t m;
    for (m = 0; m < M; ++m) {
      A[row][m] *= unit;
    }
  }
}

/* row operation of type 2 will take two matrices, A and P,
and transform them the same way.
For parameters: point_ters to both matrices,
the row to add to, the row to add from,
and the multiple of the row we are adding
*/

// adds mult*row addFrom to row addTo for column size M
// transpose of this operation is flipping addTo and addFrom
// inverse of this operation is changing mult to -1*mult
void type2rowM(std::vector<std::vector<int_t>> &A, int_t addTo, int_t addFrom,
               int_t mult) {
  int_t M = A.size();
  int_t m;
  for (m = 0; m < M; ++m) {
    A[addTo][m] += mult * A[addFrom][m];
  }
}

// int_terchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowM(std::vector<std::vector<int_t>> &A, int_t row1, int_t row2) {
  int_t M = A.size();
  int_t m;
  for (m = 0; m < M; ++m) {
    // this swaps entries without creating a temp
    A[row1][m] += A[row2][m];
    A[row2][m] = A[row1][m] - A[row2][m];
    A[row1][m] -= A[row2][m];
  }
}

// transpose of this operation is itself
// inverse of this operation is itself
void type1rowN(std::vector<std::vector<int_t>> &A, int_t row, int_t unit) {
  int_t N = A.size();
  if (unit != 1 && unit != -1) {
    printf(
        "A type 1 row operation did not execute because the input was not a "
        "unit.");
    return;
  } else {
    int_t n;
    for (n = 0; n < N; ++n) {
      A[row][n] *= unit;
    }
  }
}

// adds mult*row addFrom to row addTo for column size M
// transpose of this operation is flipping addTo and addFrom
// inverse of this operation is changing mult to -1*mult
void type2rowN(std::vector<std::vector<int_t>> &A, int_t addTo, int_t addFrom,
               int_t mult) {
  int_t N = A.size();
  int_t n;
  for (n = 0; n < N; ++n) {
    A[addTo][n] += mult * A[addFrom][n];
  }
}

// int_terchanges row1 and row2 for column size M
// transpose of this operation is itself
// inverse of this operation is itself
void type3rowN(std::vector<std::vector<int_t>> &A, int_t row1, int_t row2) {
  int_t N = A.size();
  int_t n;
  for (n = 0; n < N; ++n) {
    // this swaps entries without creating a temp
    A[row1][n] += A[row2][n];
    A[row2][n] = A[row1][n] - A[row2][n];
    A[row1][n] -= A[row2][n];
  }
}

// column operations

// multiplies entire row by unit
void type1col(std::vector<std::vector<int_t>> &A, int_t len, int_t col,
              int_t unit) {
  if (unit != 1 && unit != -1) {
    printf(
        "A type 1 column operation did not execute because the input was not a "
        "unit.");
    return;
  } else {
    int_t i;
    for (i = 0; i < len; ++i) {
      A[i][col] *= unit;
    }
  }
}

// adds mult*column addFrom to column addTo
void type2col(std::vector<std::vector<int_t>> &A, int_t len, int_t addTo,
              int_t addFrom, int_t mult) {
  int_t i;
  for (i = 0; i < len; ++i) {
    A[i][addTo] += mult * A[i][addFrom];
  }
}

// int_terchanges col1 and col2
void type3col(std::vector<std::vector<int_t>> &A, int_t len, int_t col1,
              int_t col2) {
  int_t i;
  for (i = 0; i < len; ++i) {
    // this swaps entries without creating a temp
    A[i][col1] += A[i][col2];
    A[i][col2] = A[i][col1] - A[i][col2];
    A[i][col1] -= A[i][col2];
  }
}

// takes transpose of a size N*N matrix
void transposeN(std::vector<std::vector<int_t>> &A) {
  int_t N = A.size();
  int_t i, j;
  for (i = 0; i < N; ++i) {
    for (j = i + 1; j < N; ++j) {
      // this swaps entries without creating a temp
      A[i][j] += A[j][i];
      A[j][i] = A[i][j] - A[j][i];
      A[i][j] -= A[j][i];
    }
  }
}

void transposeM(std::vector<std::vector<int_t>> &A) {
  int_t M = A.size();
  int_t i, j;
  for (i = 0; i < M; ++i) {
    for (j = i + 1; j < M; ++j) {
      // this swaps entries without creating a temp
      A[i][j] += A[j][i];
      A[j][i] = A[i][j] - A[j][i];
      A[i][j] -= A[j][i];
    }
  }
}

// returns 1 if x is in A, 0 if x is not in A
int_t contains(std::vector<int_t> A, int_t len, int_t x) {
  int_t i;
  int_t ret = 0;
  for (i = 0; i < len; ++i) {
    if (A[i] == x) {
      ret = 1;
      break;
    }
  }
  return ret;
}

// compares objects for sorting
int_t comp(const void *pa, const void *pb) {
  int_t a = *(const int_t *)pa;
  int_t b = *(const int_t *)pb;
  if (a == b) {
    return 0;
  } else if (a < b) {
    return -1;
  } else {
    return 1;
  }
}

int_t done(std::vector<int_t> A, std::vector<int_t> B, int_t lenA, int_t lenB) {
  int_t i;
  int_t ret = 1;
  for (i = 0; i < lenA; ++i) {
    if (contains(A, lenA, i) == 0) {
      ret = 0;
      break;
    }
  }
  if (ret == 1) {
    for (i = 0; i < lenB; ++i) {
      if (contains(B, lenB, i) == 0) {
        ret = 0;
        break;
      }
    }
  }
  return ret;
}

void updateFinishedRows(std::vector<std::vector<int_t>> &A,
                        std::vector<int_t> &finishedRows) {
  int_t N = A.size();
  int_t M = N;
  int_t tempCount, m, n;
  for (n = 0; n < N; ++n) {
    tempCount = 0;
    for (m = 0; m < M; ++m) {
      if (A[n][m] != 0) {
        ++tempCount;
      }
    }
    if (tempCount < 2) {
      finishedRows[n] = n;
    } else {
      finishedRows[n] = -1;
    }
  }
}

void updateFinishedColumns(std::vector<std::vector<int_t>> &A,
                           std::vector<int_t> &finishedColumns) {
  int_t N = A.size();
  int_t M = N;
  int_t tempCount, m, n;
  for (m = 0; m < M; ++m) {
    tempCount = 0;
    for (n = 0; n < N; ++n) {
      if (A[n][m] != 0) {
        ++tempCount;
      }
    }
    if (tempCount < 2) {
      finishedColumns[m] = m;
    } else {
      finishedColumns[m] = -1;
    }
  }
}

void findLeastEntry(std::vector<std::vector<int_t>> &A,
                    std::vector<int_t> &finishedRows,
                    std::vector<int_t> &finishedColumns, int_t *tempN,
                    int_t *tempM, int_t *finished) {
  int_t N = A.size();
  int_t M = N;
  int_t m, n;
  int_t boole = 0;
  *finished = 0;
  int_t tempMin = -1;

  for (n = 0; n < N; ++n) {
    for (m = 0; m < M; ++m) {
      if (!(contains(finishedRows, N, n) == 1 &&
            contains(finishedColumns, M, m) == 1)) {
        if (A[n][m] != 0) {
          *tempN = n;
          *tempM = m;
          tempMin = absolu(A[n][m]);
          boole = 1;
          break;
        }
        if (boole == 1) {
          break;
        }
      }
    }
  }
  if (tempMin == -1) {
    printf("There is no valid least entry\n");
    *finished = 1;
    return;
  }
  for (n = 0; n < N; ++n) {
    for (m = 0; m < M; ++m) {
      if (contains(finishedRows, N, n) == 0 ||
          contains(finishedColumns, M, m) == 0) {
        if (A[n][m] != 0) {
          if (absolu(A[n][m]) < tempMin) {
            *tempN = n;
            *tempM = m;
            tempMin = absolu(A[n][m]);
          }
        }
      }
    }
  }
  if (*tempM == -1 || *tempN == -1) {
    *finished = 1;
  }
}

// returns 1 if the entry divides its entire row and col; 0 otherwise
int_t dividesRowAndCol(std::vector<std::vector<int_t>> &A, int_t tempN,
                       int_t tempM) {
  int_t N = A.size();
  int_t M = N;
  int_t m, n;
  int_t ret = 1;
  if (A[tempN][tempM] == 0) {
    printf("You have called this with A[i][j] = 0\n");
    return 0;
  }
  if (tempN == -1 || tempM == -1) {
    printf("You have called this with tempN = -1 or tempM = -1\n");
    return 0;
  }
  for (n = 0; n < N; ++n) {
    if (n != tempN) {
      if (A[n][tempM] % A[tempN][tempM] != 0) {
        ret = 0;
      }
    }
  }

  for (m = 0; m < M; ++m) {
    if (m != tempM) {
      if (A[tempN][m] % A[tempN][tempM] != 0) {
        ret = 0;
      }
    }
  }
  return ret;
}

// returns the 'q' in a = bq + r
int_t eucDiv(int_t a, int_t b) {
  if (a == 0) {
    return 0;
  }
  if (b == 0) {
    printf("You are calling a division by 0\n");
    return a;
  }
  int_t temp = a % b;
  if (temp < 0) {
    temp += absolu(b);
  }
  return (a - temp) / b;
}

void findLeastEntry2(std::vector<std::vector<int_t>> &A, int_t *tempN,
                     int_t *tempM, int_t counter) {
  int_t N = A.size();
  int_t M = N;
  int_t m, n;
  int_t tempMin = -1;
  int_t boole = 0;
  if (counter >= rank_) {
    return;
  } else {
    for (n = counter; n < N; ++n) {
      for (m = counter; m < M; ++m) {
        if (A[n][m] != 0) {
          *tempN = n;
          *tempM = m;
          tempMin = absolu(A[n][m]);
          boole = 1;
          break;
        }
        if (boole == 1) {
          break;
        }
      }
    }
    if (tempMin == -1) {
      printf("This should not execute Error in findLeastEntry2.\n");
      return;
    }

    for (n = counter; n < N; ++n) {
      for (m = counter; m < M; ++m) {
        if (A[n][m] != 0) {
          if (A[n][m] < tempMin) {
            *tempN = n;
            *tempM = m;
            tempMin = A[n][m];
          }
        }
      }
    }

    if (*tempN == -1 || *tempM == -1) {
      printf("tempN or tempN is equal to -1\n");
    }
  }
}

// combines all type 1 row operations int_to one function
void rowOperations1(std::vector<std::vector<int_t>> &A,
                    std::vector<std::vector<int_t>> &P,
                    std::vector<std::vector<int_t>> &Pinv, int_t n,
                    int_t unit) {
  type1rowM(A, n, unit);
  type1rowN(P, n, unit);
  type1rowN(Pinv, n, unit);
}

// combines all type 2 row operations int_to one function
void rowOperations2(std::vector<std::vector<int_t>> &A,
                    std::vector<std::vector<int_t>> &P,
                    std::vector<std::vector<int_t>> &Pinv, int_t n, int_t tempN,
                    int_t tempMult) {
  type2rowM(A, n, tempN, -1 * tempMult);
  type2rowN(P, n, tempN, -1 * tempMult);
  type2rowN(Pinv, tempN, n, tempMult);
}

// combines all type 3 row operations int_to one function
void rowOperations3(std::vector<std::vector<int_t>> &A,
                    std::vector<std::vector<int_t>> &P,
                    std::vector<std::vector<int_t>> &Pinv, int_t row1,
                    int_t row2) {
  if (row1 != row2) {
    // row2 < min(N, M) && row1 < min(N, M) &&
    type3rowM(A, row1, row2);
    type3rowN(P, row1, row2);
    type3rowN(Pinv, row1, row2);
  }
}

// combines all type 2 column operations int_to one function
void columnOperations2(std::vector<std::vector<int_t>> &A,
                       std::vector<std::vector<int_t>> &Q,
                       std::vector<std::vector<int_t>> &Qinv, int_t m,
                       int_t tempM, int_t tempMult) {
  int_t N = A.size();
  int_t M = N;
  type2col(A, N, m, tempM, -1 * tempMult);
  type2col(Q, M, m, tempM, -1 * tempMult);
  type2col(Qinv, M, tempM, m, tempMult);
}

// combines all type 2 column operations int_to one function
void columnOperations3(std::vector<std::vector<int_t>> &A,
                       std::vector<std::vector<int_t>> &Q,
                       std::vector<std::vector<int_t>> &Qinv, int_t col1,
                       int_t col2) {
  int_t N = A.size();
  int_t M = N;
  if (col1 != col2) {
    type3col(A, N, col1, col2);
    type3col(Q, M, col1, col2);
    type3col(Qinv, M, col1, col2);
  }
}

void smithTransform(std::vector<std::vector<int_t>> &A,
                    std::vector<std::vector<int_t>> &P,
                    std::vector<std::vector<int_t>> &Pinv,
                    std::vector<std::vector<int_t>> &Q,
                    std::vector<std::vector<int_t>> &Qinv, int_t pos) {
  int_t tempMult;
  columnOperations2(A, Q, Qinv, pos, pos + 1, -1);

  while (A[pos + 1][pos] != 0) {
    tempMult = eucDiv(A[pos + 1][pos], A[pos][pos]);

    // runs Euclidean division
    rowOperations2(A, P, Pinv, pos + 1, pos, tempMult);

    // makes A[pos][pos] entry larger than A[pos][pos+1]
    if (A[pos + 1][pos] != 0) {
      rowOperations3(A, P, Pinv, pos, pos + 1);
    }
  }
  tempMult = eucDiv(A[pos][pos + 1], A[pos][pos]);
  columnOperations2(A, Q, Qinv, pos + 1, pos, tempMult);
  if (A[pos][pos + 1] != 0) {
    printf("smithTransform failed. This should never happen");
  }
}

void makeAllDiagsPositive(std::vector<std::vector<int_t>> &A,
                          std::vector<std::vector<int_t>> &P,
                          std::vector<std::vector<int_t>> &Pinv) {
  int_t N = A.size();
  int_t M = N;
  int_t n, m;
  for (n = 0; n < N; ++n) {
    for (m = 0; m < M; ++m) {
      if (A[n][m] < 0) {
        rowOperations1(A, P, Pinv, n, -1);
      }
    }
  }
}

void orderDiagonals(std::vector<std::vector<int_t>> &A,
                    std::vector<std::vector<int_t>> &P,
                    std::vector<std::vector<int_t>> &Pinv,
                    std::vector<std::vector<int_t>> &Q,
                    std::vector<std::vector<int_t>> &Qinv) {
  int_t counter = 0;
  int_t tempM, tempN;
  while (counter < rank_) {
    tempM = -1;
    tempN = -1;
    findLeastEntry2(A, &tempN, &tempM, counter);
    if (tempN != counter || tempM != counter) {
      rowOperations3(A, P, Pinv, tempN, counter);
      columnOperations3(A, Q, Qinv, tempM, counter);
    }
    ++counter;
  }
}

// returns first n where A[n+1][n+1] % A[n][n] != 0. 0 if in smith form
int_t checkSmith(std::vector<std::vector<int_t>> &A) {
  int_t smith = -1;
  int_t n;
  for (n = 0; n < rank_ - 1; ++n) {
    if (A[n + 1][n + 1] % A[n][n] != 0) {
      smith = n;
      break;
    }
  }
  return smith;
}

// calculates the rank_ of a diagonal matrix
int_t getRank(std::vector<std::vector<int_t>> &A) {
  int_t N = A.size();
  int_t M = N;
  int_t rank_ = 0;
  int_t n, m;
  for (n = 0; n < N; ++n) {
    for (m = 0; m < M; ++m) {
      if (A[n][m] != 0) {
        ++rank_;
      }
    }
  }
  return rank_;
}

// slight customization for Prof. Zhou's prime factorization algorithm
void primeFactorize(int_t Num, std::vector<int_t> &prime) {
  int_t i, j, k;
  int_t index;
  int_t start;
  std::vector<int_t> primePower(maxED);
  initializeZero(primePower, maxED);
  int_t tempEntry;
  for (j = 0; j < maxED; ++j) {
    if (prime[j] == 0) {
      index = j;
      start = j;
      break;
    }
  }

  // check if 2 is a factor
  if (Num > 1 && Num % 2 == 0) {
    prime[index] = 2;
  }

  // take out all factors of 2
  while (Num > 1 && Num % 2 == 0) {
    primePower[index] += 1;
    Num = Num / 2;
  }

  while (Num > 1 && i > 0 && 4 * i * i < Num) {
    if (Num % (2 * i + 1) == 0) {
      index++;
      prime[index] = 2 * i + 1;
    }
    while (Num > 1 && Num % (2 * i + 1) == 0) {
      primePower[index]++;
      Num = Num / (2 * i + 1);
    }
    ++i;
  }

  // it should actually never hit here...
  if (Num > 1) {
    index++;
    primePower[index] = 1;
    prime[index] = Num;
  }

  for (j = start; j < index + 1; ++j) {
    tempEntry = prime[j];
    for (k = 1; k < primePower[j]; ++k) {
      prime[j] *= tempEntry;
    }
  }
}