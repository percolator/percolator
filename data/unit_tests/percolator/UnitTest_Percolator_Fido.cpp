/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
/*
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/* This file include test cases for the EludeCaller class */
#include <gtest/gtest.h>

#include "PackedVector.cpp"
#include "PackedMatrix.cpp"
#include "Set.cpp"
#include "Numerical.cpp"
#include "Vector.cpp"
#include "BaseSpline.h"

class FidoVectorTest : public ::testing::Test {
 protected:
   virtual void SetUp() {
     pv = PackedVector();
     pv.packedAddElement(1,10);
     pv.packedAddElement(3,20);
     pv.packedAddElement(5,30);
     pv2 = PackedVector();
     pv2.packedAddElement(2,5);
     pv2.packedAddElement(5,30);
     pv2.packedAddElement(9,-10);
   }
   virtual void TearDown() {}

   PackedVector pv;
   PackedVector pv2;
};

TEST_F(FidoVectorTest, PackedVector){
  pv = PackedVector();
  EXPECT_EQ(0,pv.numberEntries());
  EXPECT_EQ(0,pv.size());

  pv = PackedVector(3);
  EXPECT_EQ(3,pv.numberEntries());
  EXPECT_EQ(3,pv.size());
  for(int i=0; i<pv.size(); i++){
    EXPECT_EQ(i,pv.index(i));
    EXPECT_EQ(0,pv[i]);
  }

  pv = PackedVector(3,1);
  EXPECT_EQ(3,pv.numberEntries());
  EXPECT_EQ(3,pv.size());
  for(int i=0; i<pv.size(); i++){
    EXPECT_EQ(i,pv.index(i));
    EXPECT_EQ(1,pv[i]);
  }
}

TEST_F(FidoVectorTest, packedReplace){
  pv = PackedVector(3,1);
  int size = pv.size();
  int numberEntries = pv.numberEntries();

  int ind = 1;
  double value = 5;
  pv.packedReplace(ind,value);
  EXPECT_EQ(value,pv[ind]);
  EXPECT_EQ(size,pv.size());
  EXPECT_EQ(numberEntries,pv.numberEntries());
}

TEST_F(FidoVectorTest, packedAddElement){
  pv = PackedVector(3,1);
  int size = pv.size();
  int numberEntries = pv.numberEntries();

  int ind = 10;
  double value = 5;
  pv.packedAddElement(ind,value);
  EXPECT_EQ(numberEntries+1, pv.numberEntries());
  EXPECT_EQ(size+1, pv.size());
  EXPECT_EQ(value,pv[pv.find(ind)]);
}

TEST_F(FidoVectorTest, swapElements){
  EXPECT_EQ(10,pv[0]);
  EXPECT_EQ(20,pv[1]);
  pv.swapElements(0,1);
  EXPECT_EQ(10,pv[1]);
  EXPECT_EQ(20,pv[0]);
}

TEST_F(FidoVectorTest, makeSparse){
  EXPECT_EQ(3,pv.size());
  EXPECT_EQ(3,pv.numberEntries());
  PackedVector sparse = pv.makeSparse();
  EXPECT_EQ(6,sparse.size());
  EXPECT_EQ(6,sparse.numberEntries());
}

TEST_F(FidoVectorTest, packedProd){
  const PackedVector old = pv;
  pv.packedProd(-1);
  EXPECT_EQ(-10,pv[pv.find(1)]);
  EXPECT_EQ(-20,pv[pv.find(3)]);
  EXPECT_EQ(-30,pv[pv.find(5)]);
  pv.packedProd(-1);
  EXPECT_TRUE(pv == old);
}

TEST_F(FidoVectorTest, packedSubtract){
  const PackedVector pv_const = pv;
  const PackedVector pv2_const = pv2;
  PackedVector diff = pv_const.packedSubtract(pv2_const);
  EXPECT_EQ(10,diff[diff.find(1)]);
  EXPECT_EQ(-5,diff[diff.find(2)]);
  EXPECT_EQ(20,diff[diff.find(3)]);
  EXPECT_EQ(10,diff[diff.find(9)]);
  EXPECT_EQ(4,diff.numberEntries());
}

TEST_F(FidoVectorTest, packedAdd){
  const PackedVector pv_const = pv;
  const PackedVector pv2_const = pv2;
  PackedVector diff = pv_const.packedAdd(pv2_const);
  EXPECT_EQ(10,diff[diff.find(1)]);
  EXPECT_EQ(5,diff[diff.find(2)]);
  EXPECT_EQ(20,diff[diff.find(3)]);
  EXPECT_EQ(60,diff[diff.find(5)]);
  EXPECT_EQ(-10,diff[diff.find(9)]);
}

TEST_F(FidoVectorTest, packedDotProd){
  const PackedVector pv_const = pv;
  const PackedVector other_const = pv;
  double prod = pv_const.packedDotProd(other_const);
  EXPECT_EQ(1400,prod);
}

class FidoMatrixTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    pm = PackedMatrix(3,3);
    pm[0].packedAddElement(0,1);
    pm[0].packedAddElement(1,3);
    pm[0].packedAddElement(2,-2);
    pm[1].packedAddElement(0,3);
    pm[1].packedAddElement(1,5);
    pm[1].packedAddElement(2,6);
    pm[2].packedAddElement(0,2);
    pm[2].packedAddElement(1,4);
    pm[2].packedAddElement(2,3);
    res = PackedVector();
    res.packedAddElement(0,5);
    res.packedAddElement(1,7);
    res.packedAddElement(2,8);
   }
   virtual void TearDown() {}

   PackedMatrix pm;
   PackedMatrix pm2;
   PackedVector res;
};

TEST_F(FidoMatrixTest, PackedMatrix){
  EXPECT_EQ(3,pm.numCols());
  EXPECT_EQ(3,pm.numRows());
  EXPECT_EQ(1,pm[0][0]);
  EXPECT_EQ(3,pm[0][1]);
  EXPECT_EQ(-2,pm[0][2]);
  EXPECT_EQ(3,pm[1][0]);
  EXPECT_EQ(5,pm[1][1]);
  EXPECT_EQ(6,pm[1][2]);
  EXPECT_EQ(2,pm[2][0]);
  EXPECT_EQ(4,pm[2][1]);
  EXPECT_EQ(3,pm[2][2]);
  pm = PackedMatrix();
  EXPECT_EQ(0,pm.numCols());
  EXPECT_EQ(0,pm.numRows());
  SetUp();
}

TEST_F(FidoMatrixTest, solveInPlace){
BaseSpline::solveInPlace(pm,res);
EXPECT_EQ(-15,(int)res[0]);
EXPECT_EQ(8,(int)res[1]);
EXPECT_EQ(2,(int)res[2]);
}
