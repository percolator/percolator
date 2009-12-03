#include <LTSRegression.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <float.h>

// initialize static variables
int LTSRegression::noSubsets = 3;
// maximum difference btw q2 and q1 to achieve convergence
double LTSRegression::epsilon = 0.0001;

// function to compare 2 data points according to r
bool compareDataPoints (dataPoint x, dataPoint y) { return (x.absr < y.absr); }

LTSRegression::LTSRegression(): h(-1), regCoefficients(make_pair(0.0, 0.0))
{
}

LTSRegression::~LTSRegression()
{
}

void LTSRegression::setData(vector<double> & x, vector<double> & y)
{
	dataPoint tmp;
	for(int i = 0; i < x.size(); ++i)
	{
		tmp.x = x[i];
		tmp.y = y[i];
		tmp.absr = -1.0;
		data.push_back(tmp);
	}

	// since we definitely expect less than 25% contamination, we set h as 0.75*n
	h = (int) round(0.75 * x.size());
}

// for our case, constructing a random p-subset is equivalent to build the equation of a line through 2 randomly
// chosen points
vector<dataPoint> LTSRegression::getInitialHSubset()
{
	double a, b;
	int i1, i2;
	int n = data.size();

	// generate two random indices
	i1 = rand() % n;
	i2 = rand() % n;
	while (i2 == i1)
		i2 = rand() % n;

	// calculate the a and b of the equation of the line going through the two points selected above (y = ax + b)
	a = (data[i2].y - data[i1].y) / (data[i2].x - data[i1].x);
	b = data[i1].y - (data[i1].x * a);

	// calculate the residuals and sort the data according to these values
	fillResiduals(make_pair(a,b));
	partial_sort(data.begin(), data.begin() + h, data.end(), compareDataPoints);

	return vector<dataPoint>(data.begin(), data.begin() + h);
}

// fill the absolute values of the residuals
void LTSRegression::fillResiduals(pair<double, double> par)
{
	for (int i = 0; i < data.size(); ++i)
		data[i].absr = abs(data[i].y - (par.first * data[i].x + par.second));
}

// fit a line to the points in h using the least-squares method
pair<double, double> LTSRegression::fitLSLine(vector<dataPoint> hdata)
{
	double sumxy = 0.0, sumx = 0.0, sumy = 0.0, sumxsq = 0.0;
	double a, b;

	for(int i = 0; i < h; ++i)
	{
		sumx += hdata[i].x;
		sumy += hdata[i].y;
		sumxy += hdata[i].x * hdata[i].y;
		sumxsq += pow(hdata[i].x,2);
	}

	a = ((h * sumxy) - (sumx * sumy)) / ((h * sumxsq) - (pow(sumx, 2)));
	b = (sumy / h) - (a * (sumx / h));

	//cout << "a , b = " << a << ", " << b <<  endl;
	return make_pair(a,b);
}

vector<dataPoint> LTSRegression::performCstep(vector<dataPoint> hOld)
{
	pair<double, double> par;

	// fit a least-squares regression line using data points from hold
	par = fitLSLine(hOld);
	// calculate and fill the residuals for all the data points
	fillResiduals(par);
	// sort the data points according to absolute values of residuals
	partial_sort(data.begin(), data.begin() + h, data.end(), compareDataPoints);

	return vector<dataPoint> (data.begin(), data.begin() + h);
}

// calculate the sum of squared residuals
double LTSRegression::calculateQ()
{
	double res = 0.0;

	for(int i = 0; i < h; ++i)
		res += pow(data[i].absr, 2);

	return res;
}

// get the 3 subsets that have the lowest Q (sum of squared residuals)
vector< vector<dataPoint> > LTSRegression::getBestHSubsets()
{
	vector< vector<dataPoint> > res;
	partial_sort(data.begin(), data.begin() + h + 2, data.end(), compareDataPoints);
	vector<dataPoint> v1(data.begin(), data.begin() + h);
	vector<dataPoint> v2(data.begin(), data.begin() + h - 1);
	vector<dataPoint> v3(data.begin(), data.begin() + h - 2);
	v2.push_back(data[h]);

	if ((pow(data[h + 1].absr,2) - pow(data[h-1].absr,2)) <= (pow(data[h].absr,2) - pow(data[h-2].absr,2)))
	{
		v3.push_back(data[h-2]);
		v3.push_back(data[h+1]);
	}
	else
	{
		v3.push_back(data[h]);
		v3.push_back(data[h-1]);
	}

	res.push_back(v1);
	res.push_back(v2);
	res.push_back(v3);

	return res;
}

void LTSRegression::runLTS()
{
	double q1, q2, bestq = DBL_MAX;
	pair<double, double> par;
	vector<dataPoint> hold, hnew, besth(data.begin(), data.end() + h);
	vector <vector<dataPoint> > bestSubsets;

	srand ( time(NULL) );

	for(int i = 0; i < noSubsets; ++i)
	{
		// get the first subset (data will be sorted according to abs(residuals)
		hold = getInitialHSubset();

		// apply 2 C-steps
		hnew = performCstep(hold);
		hold = performCstep(hnew);

		// take the first 3 subsets with lowest Q
		bestSubsets = getBestHSubsets();

		// for each such subset perform C-steps until convergence
		for(int j = 0; j < bestSubsets.size(); ++j)
		{
			hold = bestSubsets[j];
			par = fitLSLine(hold);
			fillResiduals(par);
			q2 = calculateQ();
			do
			{
				q1 = q2;
				hnew = performCstep(hold);
				par = fitLSLine(hnew);
				fillResiduals(par);
				q2 = calculateQ();
				hold = hnew;
			}
			while( abs(q2 - q1) > epsilon);

			if (q2 < bestq)
			{
				regCoefficients = par;
				bestq = q2;
				besth = hnew;
			}
		}
	}
	printDataPoints();
	cout << "Final equation: y = " << regCoefficients.first << " * x + " << regCoefficients.second << endl;
	cout << "Data points used to generate this equation: " << endl;
	printVector(besth);

	cout << "--------------------------------------\n";
	vector<double> d;
	for(int i = 0; i < data.size(); ++i)
		d.push_back(data[i].x);
	vector<double> preds;
	preds = predict(d);
	for(int i = 0; i < data.size(); ++i)
		cout << data[i].x << " " << data[i].y << "    " << preds[i] << endl;
}

vector<double> LTSRegression::predict(vector<double> & x)
{
	vector<double> predictions;
	double a = regCoefficients.first, b = regCoefficients.second;

	for(int i = 0; i < x.size(); ++i)
		predictions.push_back((a * x[i]) + b);

	return predictions;
}

void LTSRegression::printVector(vector<dataPoint> v)
{
	for(int i = 0; i < v.size(); ++i)
		cout << v[i].x << " " << v[i].y << " " << v[i].absr << endl;
	cout << endl;
}


void LTSRegression::printDataPoints()
{
	for(int i = 0; i < h; ++i)
		cout << data[i].x << " " << data[i].y << " " << data[i].absr << endl;

	cout << "----" << endl;

	for(int i = h; i < data.size(); ++i)
		cout << data[i].x << " " << data[i].y << " " << data[i].absr << endl;
}
