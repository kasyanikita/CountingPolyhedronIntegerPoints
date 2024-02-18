#ifndef COUNTINGINTEGERPOINTS_SNF_H_
#define COUNTINGINTEGERPOINTS_SNF_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "global_defs.h"

using namespace GroupIP;

void SNF(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &);

void type1rowM(std::vector<std::vector<int_t>> &, int_t, int_t);
void type2rowM(std::vector<std::vector<int_t>> &, int_t, int_t, int_t);
void type3rowM(std::vector<std::vector<int_t>> &, int_t, int_t);

void type1rowN(std::vector<std::vector<int_t>> &, int_t, int_t);
void type2rowN(std::vector<std::vector<int_t>> &, int_t, int_t, int_t);
void type3rowN(std::vector<std::vector<int_t>> &, int_t, int_t);

void type1col(std::vector<std::vector<int_t>> &, int_t, int_t, int_t);
void type2col(std::vector<std::vector<int_t>> &, int_t, int_t, int_t, int_t);
void type3col(std::vector<std::vector<int_t>> &, int_t, int_t, int_t);

void rowOperations1(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, int_t, int_t);

void columnOperations2(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					   std::vector<std::vector<int_t>> &, int_t, int_t, int_t);
void rowOperations2(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, int_t, int_t, int_t);

void rowOperations3(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, int_t, int_t);
void columnOperations3(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					   std::vector<std::vector<int_t>> &, int_t, int_t);

void transposeN(std::vector<std::vector<int_t>> &);
void transposeM(std::vector<std::vector<int_t>> &);

int_t contains(std::vector<int_t>, int_t, int_t);
int_t done(std::vector<int_t>, std::vector<int_t>, int_t, int_t);
void makeAllDiagsPositive(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
						  std::vector<std::vector<int_t>> &);

void updateFinishedRows(std::vector<std::vector<int_t>> &, std::vector<int_t> &);
void updateFinishedColumns(std::vector<std::vector<int_t>> &, std::vector<int_t> &);

int_t dividesRowAndCol(std::vector<std::vector<int_t>> &, int_t, int_t);
int_t eucDiv(int_t, int_t);
int_t min(int_t, int_t);
int_t comp(const void *a, const void *b);

int_t checkSmith(std::vector<std::vector<int_t>> &);
int_t getRank(std::vector<std::vector<int_t>> &);

void smithTransform(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, int_t);
void orderDiagonals(std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &, std::vector<std::vector<int_t>> &,
					std::vector<std::vector<int_t>> &);

void findLeastEntry(std::vector<std::vector<int_t>> &, std::vector<int_t> &,
					std::vector<int_t> &, int_t *, int_t *, int_t *);
void findLeastEntry2(std::vector<std::vector<int_t>> &, int_t *, int_t *, int_t);

void primeFactorize(int_t, std::vector<int_t> &);
void initializeZero(std::vector<int_t> &, int_t);

#endif  // COUNTINGINTEGERPOINTS_TODDFFT_H_