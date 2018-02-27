#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <math.h>     
using namespace std;
 
struct Point {
    vector<double> coordinates;
};
 
typedef std::vector<Point> DataSet;
typedef DataSet::iterator DataSet_iterator;
 
double dot_product(Point a, Point b){ // computing the dot product of two vectors of doubles
	int dimension=a.coordinates.size();	
	double result=0.;
	for (int k=0; k<dimension; k++){
		result += a.coordinates[k]*b.coordinates[k];
	}
	return result;
} 

struct compare_along_direction {
	compare_along_direction(Point direction){this->direction=direction;};
	bool operator ()(const Point& left, const Point& right) const {
		return dot_product(left, direction) < dot_product(right, direction);
	}
	Point direction;
};

bool compare_along_direction_simple(Point left, Point right, Point direction){
	return dot_product(left, direction) < dot_product(right, direction);
}

// returns a distance matrix
vector<vector<double>> compute_all_distances(DataSet data){
	
	int number_of_points=data.size();
	int dimension=data[0].coordinates.size();
	
	vector<vector<double>> distance_matrix (number_of_points, vector<double>(number_of_points, 0));
	
	for (int i=0; i < number_of_points; i++) {
		
		for (int j=0; j < i; j++){
			
			vector<double> point_i=data[i].coordinates;
			vector<double> point_j=data[j].coordinates;
			
			double distance=0.;
			
			for(int k=0; k < dimension; k++){
				double t = point_j[k] - point_i[k];
				distance += t*t;
			}
			
			distance_matrix[i][j]=distance;
			distance_matrix[j][i]=distance;
		}
	}
	return distance_matrix;
}
// finding the maximum distance within a distance matrix
double find_maximum_distance(vector<vector<double>> distance_matrix){
	
	double maximum=0.;
	int number_of_points=distance_matrix.size();
	
	for(int i=0; i<number_of_points; i++){
		for(int j=0; j<i; j ++){
			if(distance_matrix[i][j] > maximum){
				maximum = distance_matrix[i][j];
			}
		}
	}
	return sqrt(maximum);
}
// wrapping both distance computation and maximum exploration
double compute_diameter(DataSet data){
	vector<vector<double>> distance_matrix=compute_all_distances(data);
	double diameter=find_maximum_distance(distance_matrix);
	return diameter;
}

//generating the canonical base
DataSet canonical_base(int dimension){
	DataSet base;
	Point zeros;
	for(int k=0; k<dimension; k++){
		zeros.coordinates.push_back(0.);
	}
	for(int i = 0; i<dimension; i++){
		zeros.coordinates[i] = 1.;
		base.push_back(zeros);
		zeros.coordinates[i] = 0.;
	}
	
	return base;
}

 
struct Node {
    Point point;
    Node *left;
    Node *right;
    double diameter;
    Node(DataSet dataset, DataSet base, int index_direction, int depth);
};
 
Node::Node( DataSet dataset_input, DataSet base, int index_direction, int depth = 0) {
	// defining the beginning and ending iterators
	DataSet_iterator begin = dataset_input.begin();
	DataSet_iterator end = dataset_input.end();
	
	depth = depth;
	// terminal cases
    left = right = 0;
	    if (begin == end){
        return ;
    }
    if( end - begin == 1) {
        point = *begin;
        return;
    }
    
    // get the direction for the next split
    Point new_direction = base[index_direction];
    // get the dimension of the points
    int dimension = new_direction.coordinates.size();
    // sort the dataset along the chosen direction
    std::sort( begin, end, compare_along_direction(new_direction));
    // compute the diameter of the dataset
    diameter = compute_diameter(dataset_input);
    
    // computing the median
    DataSet_iterator median = begin + (end - begin)/2;
    Point median_point;
    median_point = *(median);
    
    // computing the value of the jittered split along the chosen direction
    double random_number = (double)rand()/RAND_MAX;
    random_number = -1. + random_number*2.;
    cout << "random number : " << random_number << endl;
    double jittered_split_value = diameter*3. / sqrt((double)dimension);
    
    jittered_split_value += dot_product(new_direction, median_point);
    cout << "diameter " << diameter << " at depth " << depth << endl;
    cout << "split value " << jittered_split_value << endl;
    /*
    // computing an iterator to get the data into two groups
    Point jittered_point = new_direction; // this point is not in our data set but we will take the closest to it
    for(int k = 0; k<dimension; k++){ // we need to multiply the unit vector by the jittered split value
		jittered_point.coordinates[k] *= jittered_split_value;
	}
    DataSet_iterator jittered_split_point;
    jittered_split_point = std::upper_bound(begin, end, jittered_point, compare_along_direction(new_direction));
	
	if(jittered_split_point==begin){
		jittered_split_point+=1;
	} else if(jittered_split_point==end){
		jittered_split_point-=1;
	}
	*/
    // grouping the data into the subnodes
    DataSet dataset_left, dataset_right;
    
    for(DataSet_iterator pt = begin; pt!=end; pt++){
		if(dot_product(*pt, new_direction)> jittered_split_value){
			dataset_right.push_back(*pt);
		} else {
			dataset_left.push_back(*pt);
		}
	}
	
	cout << "Size of left dataset : " << dataset_left.size() << endl;
	cout << "Size of right dataset : " << dataset_right.size() << endl;
	
    // setting the value of the point
    //point = *(jittered_split_point);
    // recursivity
    if ( dataset_left.size()>0){
	    left  = new Node(dataset_left , base, (index_direction + 1) % dimension, depth + 1);
    }
    if ( dataset_right.size()>0){
	    right = new Node(dataset_right,base, (index_direction + 1) % dimension, depth + 1);
    }
}
 
int main(){
	Point A;
	Point B;
	Point C;
	Point D;
	Point E;
	Point F;
	
	A.coordinates = {1.,2.,14.};
	B.coordinates = {1.,0.,-5.};
	C.coordinates = {1.,5.,10.};
	D.coordinates = {3.,2.,-4};
	E.coordinates = {-1.,4.,1.};
	F.coordinates = {0.,0.,10.};
	
	Point E1;
	Point E2;
	Point E3;
	E1.coordinates = {1.,0.,0.};
	E2.coordinates = {0.,1.,0.};
	E3.coordinates = {0.,0.,1.};
    
	DataSet dataset = {A,B,C,D,E,F};
	DataSet base = {E1,E2, E3};
	base = canonical_base(1000);
	
	DataSet bigdataset;
	
	ifstream input("/home/paul/Desktop/random_uniform_dataset.txt");
    
    for (int i=0; i < 1000; i++) {
		Point pt;
			
		for(int j=0; j<1000; j++){
			double number;
			input >> number;
			pt.coordinates.push_back(number);
			}
		bigdataset.push_back(pt);
	}
	
	
	Node tree = Node(bigdataset,base, 0);
	
	
	return 0;	
}
