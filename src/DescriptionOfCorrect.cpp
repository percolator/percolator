#include "DescriptionOfCorrect.h"

DescriptionOfCorrect::DescriptionOfCorrect()
{
}

DescriptionOfCorrect::~DescriptionOfCorrect()
{
}

float DescriptionOfCorrect::krokhin_index['Z'-'A'+1] = 
         {0.8, 0.0, -0.8, -0.5, 0.0,  10.5, -0.9, -1.3, 8.4, 0.0, -1.9,9.6,5.8,
          -1.2,0.0,0.2,-0.9,-1.3,-0.8,0.4,0.0,5.0,11.0,0.0,4.0,0.0};
float DescriptionOfCorrect::hessa_index['Z'-'A'+1] = 
         {0.11,0.0,-0.13,3.49,2.68,-0.32,0.74,2.06,-0.60,0.0,2.71,-0.55,-0.10,2.05,
          0.0,2.23,2.36,2.58,0.84,0.52,0.0,-0.31,0.30,0.0,0.68,0.0};
float DescriptionOfCorrect::kytedoolittle_index['Z'-'A'+1] =
         {1.80,0.0,2.50,-3.50,-3.50,2.80,-0.40,-3.20,0.0,4.50,3.90,3.80,1.90,-3.50,
          0.0,-1.60,-3.50,-4.50,-0.80,-0.70,0.0,4.20,-0.90,0.0,-1.30,0.0};


