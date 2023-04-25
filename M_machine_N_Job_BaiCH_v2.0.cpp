#include <iostream>
#include <fstream>    // read the file 
#include <iomanip>    // setw()
#include <vector>
#include <numeric>    // accumulate(), to get the sum of vector
#include <algorithm>  // swap(), swap 2 variability
using namespace std;
using std::vector;
using std::cout;
/*________________Define the file's name by users__________________*/   

const char* file_name = "Route_3_3_2_0.txt";
//const char* file_name = "Example_in_Paper.txt";     // Test file : Example in Paper(4 Jobs, 5 machines)
/*________________________Define Over______________________________*/

double Larger(double x1, double x2)   // This microfunction is used for MAKESPAN CALCULATING
{
	return ((x1 >= x2) ? x1 : x2);
}

/*****************This Fuction is to get the Makespan*************************/
double Makespan( vector<vector<double>> select)    
{
	vector<vector<double>> makespan = select;
	makespan[0][0] = select[0][0];
	for (decltype(select.size()) i = 1; i < select[0].size(); ++i) {         // Get the row 1 makespan
		makespan[0][i] = makespan[0][i - 1] + select[0][i];
	}
	for (decltype(select.size()) j = 1; j < select.size(); ++j) {            // Get the column 1 makespan
		makespan[j][0] = makespan[j - 1][0] + select[j][0];
	}
	for (decltype(select.size()) i = 1; i < select.size(); ++i) {            // Get another elements
		for (decltype(select.size()) j = 1; j < select[0].size(); ++j) {
			makespan[i][j] = Larger(makespan[i - 1][j], makespan[i][j - 1]) + select[i][j];
		}
	}
	double rslt = makespan[makespan.size() - 1][makespan[0].size() - 1];      // Total makespan is the last element
	return rslt;
}

/*_______________________MAIN FUNCTION_______________________*/

int main()
{
	// Define a vector in 2-Demension, Job_i to Machine_j, t_ij
	int rows, cols;
	vector<vector<double>> process_time;      // 2-Demension Matrix, t_ij = Job i in Machine j's process time
	cout << "Opening file : " << file_name << endl << endl;
	ifstream inFile(file_name);

	// Get the number of rows and number of columns
	if (inFile.is_open()) {
		inFile >> rows >> cols;
		process_time.resize(rows, vector<double>(cols, 0));
	}

	if (inFile.is_open()) {
		for (auto i = 0; i < rows; ++i) {
			for (auto j = 0; j < cols; ++j) {
				inFile >> process_time[i][j];
			}
		}
	}

	// Close the file
	inFile.close();

	// __________________________Difine Over___________________________
	std::cout << "t_ij, job_i in machine_j's process time:" << endl;
	for (vector<double> i : process_time) {
		for (double j : i) {
			cout << setw(5) << j ;
		}
		cout << endl;
	}
	cout << endl;

	/******Step1: Get the sum of Job_1 to Job_n******/

	vector<double> sum_time;       // Job1 to Jobn , Get the total time
	for (vector<double> i : process_time) {
		sum_time.push_back(
			std::accumulate(i.begin(), i.end(), 0.0)
		);
	}
	/*  This code is for testing, make sure the elements of sum_time  */
	cout << "The sums of Job_1 to Job_n are: " << endl;
	for (double sum : sum_time) {
		cout << sum << " ";
	}
	cout << endl << endl;
	
	/****** Step2: Selcet 2 larggest sum_time and corresponding Job[x] ******/

	// In this code is sort from largest to smallest
	// The vector "sequence" is an important output*********
	/*!!!!!!!!!Never Touch This!!!!!!!!*/
	vector<double> sum_time_jud;     // vector sum_time_jud is used to sequence, = sum_time
	/*!!!!!!!!!Never Touch This!!!!!!!!*/

	sum_time_jud = sum_time;

	// Vector sequence contain the sum time's index**********
    vector<decltype(sum_time_jud.size())> sequence(sum_time_jud.size());
	iota(sequence.begin(), sequence.end(), 0);     // sequence = {0,1,2......n-1}

	for (decltype(sum_time_jud.size()) i = 0; i < (sum_time_jud.size() - 1); ++i) {
		decltype(sum_time_jud.size()) max_index = i;
		for (decltype(sum_time_jud.size()) j = i + 1; j < (sum_time_jud.size()); ++j) {
			// If element[j] > max, redifine the maximum
			if (sum_time_jud[j] > sum_time_jud[max_index]) {
				max_index = j;
			}
		}
		swap(sum_time_jud[i], sum_time_jud[max_index]);
		swap(sequence[i], sequence[max_index]);
	}
	cout << "The sum time from large to small:" << endl;
	for (decltype(sum_time_jud.size()) idx : sequence) {
		cout << "process_time[" << idx << "]: " << sum_time[idx] << endl;
	}
	cout << endl;

	/********Step3 & 4*********/
	vector<vector<double>> select;
	select.push_back(process_time[sequence[0]]);

	vector<decltype(sequence.size())> select_idx;
	select_idx.push_back(sequence[0]);
	double minTime = 0;

	for (decltype(sum_time_jud.size()) i = 1; i < sequence.size(); ++i) {   // 决定插入的index和向量
		// The next vector used for inserting.(large to small)
		vector<double> newVec = process_time[sequence[i]];
		decltype(sum_time_jud.size()) newIdx = sequence[i];
		vector<vector<decltype(sum_time_jud.size())>> combination;

		// 生成插入所有可能的位置后的可能性
		// After 1st iteration : combination = {{2，0}，{0，2}}
		for (decltype(sum_time_jud.size()) j = 0; j <= select_idx.size(); ++j) {
			vector<decltype(sum_time_jud.size())> combination_idx = select_idx;
			combination_idx.insert(combination_idx.begin() + j, sequence[i]); // 生成了一个插入后的顺序
			combination.push_back(combination_idx);
		}
		vector<double> makespan_vec;    // Save the candidate makespan in combination

		/**** Iterating ****/
		cout << "Iteration "  << i << ":" << endl;
		cout << "The candidate machine sequence:" << endl;

		for (decltype(sum_time_jud.size()) a = 0; a < combination.size(); ++a) {     // 从combination[0]开始
			vector<vector<double>> time_mat;
			
			for (decltype(sum_time_jud.size()) b = 0; b < combination[0].size(); ++b) {   // 从combination[0][0] -> combination[0][max]
				cout << setw(1) <<combination[a][b] + 1 << " ";
				time_mat.push_back(process_time[combination[a][b]]);
			}
			cout << endl;
			double makespan = Makespan(time_mat);
			makespan_vec.push_back(makespan);

			
		}
		// cout the corresponding makespan in combination
		cout << "The corresponding makespan is:" << endl;
		for (double b : makespan_vec) {
			cout << b << " ";
		}
		cout << endl << endl;

		auto min_makespan = std::min_element(makespan_vec.begin(), makespan_vec.end());
		auto min_idx = std::distance(makespan_vec.begin(), min_makespan);
		minTime = *min_makespan;
		cout << "Minimum makespan: " << *min_makespan << endl
			 << "Corresponding index: " << min_idx 
			 << endl;
		select_idx = combination[min_idx];
		
		cout << endl << endl;
		
	}
	cout << "*********************************" 
		 << "*********************************"
		 << endl << endl 
	     << "The Macine sequence is : ";

	// Output the machine sequence
	for (decltype(sum_time_jud.size()) i = 0; i < (select_idx.size()-1); ++i) {
		cout << select_idx[i] + 1<< " - ";
	}
	cout << select_idx[select_idx.size() - 1]+1 << endl;

	// Output the minimum makespan
	cout << "The minimum makespan: " << minTime;

	cout << endl << endl 
		 << "*********************************" 
		 << "*********************************"
		 << endl << endl;
	
	return 0;
}