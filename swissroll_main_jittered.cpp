#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>
#include <string>
#include <math.h>
#include "matrix.h"  
using namespace std;

// structure of a point
struct Point {
    vector<double> coordinates;
};

// defining some types
typedef std::vector<Point> DataSet;
typedef DataSet::iterator DataSet_iterator; 


// structure of a node
struct Node {
    DataSet node_points;
    Node *left;
    Node *right;
    long long int node_id;
    int node_number_of_points;
    double node_diameter;
    int node_depth;
    Node(DataSet dataset, DataSet base, int index_direction, int depth, double jittered_factor, int minCellSize, long long int id);
};

// structure of the summary of a Tree
struct TreeSummary{
	vector<long long int> ids;
	vector<int> depths;
	vector<double> diameters;
	vector<int> numbers_of_points;
};

// explore a tree
void explore_tree(Node node, TreeSummary& summary){
	
	summary.ids.push_back(node.node_id);
	summary.depths.push_back(node.node_depth);
	summary.diameters.push_back(node.node_diameter);
	summary.numbers_of_points.push_back(node.node_number_of_points);
	
	if(node.left != 0){
		explore_tree(* node.left, summary);
	};
	if(node.right != 0){
		explore_tree(* node.right, summary);
	};
}
 
// dot product
double dot_product(Point a, Point b){ // computing the dot product of two vectors of doubles
	int dimension=a.coordinates.size();	
	double result=0.;
	for (int k=0; k<dimension; k++){
		result += a.coordinates[k]*b.coordinates[k];
	}
	return result;
} 

// loading data function
DataSet load_data(string path, int number_of_points, int dimension){
	DataSet output_dataset;
	ifstream input(path);
    
    for (int i=0; i < number_of_points; i++) {
		Point pt;			
		for(int j=0; j<dimension; j++){
			double number;
			input >> number;
			pt.coordinates.push_back(number);
			}
		output_dataset.push_back(pt);
	}
	
	return output_dataset;
}

// saving data function
void save_data(TreeSummary summary, string path){
	ofstream myfile (path);
	if(myfile.is_open()){
		myfile << "ID | Depth | diameter | Number of points\n" ;
		for (unsigned int i = 0; i < summary.ids.size(); i++){
			myfile << summary.ids[i];
			myfile << " | ";
			myfile << summary.depths[i];
			myfile << " | ";
			myfile << summary.diameters[i];
			myfile << " | ";
			myfile << summary.numbers_of_points[i];
			myfile << "\n";
		}
		myfile.close();
	} else {
		cout << "unable to open file" << endl;
	}
}

// compare according to a vector
struct compare_along_direction {
	compare_along_direction(Point direction){this->direction=direction;};
	bool operator ()(const Point& left, const Point& right) const {
		return dot_product(left, direction) < dot_product(right, direction);
	}
	Point direction;
};

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

// building a random orthonormal basis

// random basis 
DataSet matrix_rand_ortho(int dim){ // Generate a random orthonormal matrix with Householder transformations
    QSMatrix<double> A(dim, dim, 'i'); // The random orthnormal matrix
    QSMatrix<double> d(1,dim,0);
    srand(time(NULL));
    double number = rand()/(double)RAND_MAX;
    if (number < 0.5){d(0,dim-1)=-1;} else { d(0, dim-1)=1;}

    double norm_x=0;
    double sign=0;
    double beta = 0;
    for (int k=dim-2; k>-1; k--){
        // Generate random Householder transformation
        QSMatrix<double> x(1, dim-k, 'r');
        // norm of x
        norm_x=0;
        for (int i=0; i<dim-k; i++){
            norm_x += x(0, i)*x(0, i);
        }
        norm_x = sqrt(norm_x);

        sign = 1;
        if (x(0,0)<0){
            sign = -1;
        }
        norm_x *= sign;
        d(0,k) = -sign;
        x(0,0) += norm_x;
        beta = norm_x * x(0,0);
        // apply the transformation
        QSMatrix<double> A_(dim-k, dim, 0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                A_(i, j)=A(i+k, j);
            }
        }
        QSMatrix<double> y(1,dim,0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                y(0,j)+=x(0,i)*A_(i,j);
            }
        }
        y = y/beta;
        QSMatrix<double> outer_xy(dim-k, dim, 0);
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                outer_xy(i,j)=x(0,i)*y(0,j);
            }
        }
        A_-=outer_xy;
        for (int i=0; i<dim-k; i++){
            for (int j=0; j<dim; j++){
                A(i+k, j)=A_(i, j);
            }
        }
    }
    // change sign of rows
    for (int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            A(i,j) = A(i,j)*d(0,i);
		}
	}

    DataSet new_base;
    vector<vector<double>> values = A.mat;
    Point pt;
    
    for(int k=0; k<dim; k++){
		pt.coordinates = values[k];
		new_base.push_back(pt);
	}
    return new_base;   
}

// building the tree 
Node::Node( DataSet dataset_input, DataSet base, int index_direction, int depth = 0, double jittered_factor = 1, int minCellSize = 10, long long int id = 0) {
	node_id =id;
	node_number_of_points = dataset_input.size();
	node_depth = depth;
	// defining the beginning and ending iterators
	DataSet_iterator begin = dataset_input.begin();
	DataSet_iterator end = dataset_input.end();
	
	// terminal cases
    left = right = 0;
	    if (begin == end){
        return ;
    }
    if( end - begin <= minCellSize) {
        node_points = dataset_input;
        return;
    }
    
    // get the direction for the next split
    Point new_direction = base[index_direction];
    // get the dimension of the points
    int dimension = new_direction.coordinates.size();
    
    // sort the dataset along the chosen direction
    std::sort(begin, end, compare_along_direction(new_direction));
    // compute the diameter of the dataset
    double diameter = compute_diameter(dataset_input);
    node_diameter = diameter;
    // computing the median
    DataSet_iterator median = begin + (end - begin)/2;
    Point median_point;
    median_point = *(median);
    
    // computing the value of the jittered split along the chosen direction
    double random_number = (double)rand()/RAND_MAX;
    random_number = -1. + random_number*2.;
    //cout << "random number : " << random_number << endl;
    double jittered_split_value = diameter*3. / sqrt((double)dimension);
	jittered_split_value *=random_number * jittered_factor;
    jittered_split_value += dot_product(new_direction, median_point);
    //cout << "diameter " << diameter << " at depth " << depth << endl;
    //cout << "split value " << jittered_split_value << endl;

    // grouping the data into the subnodes
    DataSet dataset_left, dataset_right;
    
    for(DataSet_iterator pt = begin; pt!=end; pt++){
		if(dot_product(*pt, new_direction)> jittered_split_value){
			dataset_right.push_back(*pt);
		} else {
			dataset_left.push_back(*pt);
		}
	}
	
	//cout << "Size of left dataset : " << dataset_left.size() << endl;
	//cout << "Size of right dataset : " << dataset_right.size() << endl;
	
    // recursivity
    if ( dataset_left.size()>0){
	    left  = new Node(dataset_left , base, (index_direction + 1) % dimension, depth + 1, jittered_factor, minCellSize, 2* id );
    }
    if ( dataset_right.size()>0){
	    right = new Node(dataset_right,base, (index_direction + 1) % dimension, depth + 1, jittered_factor, minCellSize, 2* id + 1);
    }
}

 
 
int main(){
	// setting parameters
	int simulation_dimension = 3;
	int simulation_number_of_points = 1600;
	int simulation_number_of_iterations = 50;
	// setting output and input path
	string github_repo_path = "/home/paul/Desktop/MSc DSBA/7. Geometric Methods for Data Analysis/Github/GMDA-CS/";
	string output_path = github_repo_path + "results_swissroll_jittered/";
	string input_path = github_repo_path + "input_data/swissroll";
	// getting canonical base
	DataSet base_canonical = canonical_base(simulation_dimension);
	//getting data
	DataSet dataset = load_data(input_path, simulation_number_of_points, simulation_dimension);
	cout << "number of iterations " << simulation_number_of_iterations<< endl;
	double jittered_factor_simulation = 2;
	for(int jittered_iteration =1; jittered_iteration < 5; jittered_iteration++){
		jittered_factor_simulation *=.5;
		for(int iteration = 1; iteration <= simulation_number_of_iterations; iteration++){
			
			string iteration_string = std::to_string(iteration);
			cout << "iteration n°" << iteration << endl; 
			// defining a rotated base
			DataSet base_rotated = matrix_rand_ortho(simulation_dimension);
			// Building the two trees
			Node tree = Node(dataset,base_canonical, 0, 0, jittered_factor_simulation, 10);
			Node rotated_tree = Node(dataset,base_rotated, 0, 0, jittered_factor_simulation, 10);

			// Exploring trees
			TreeSummary tree_summary, rotated_tree_summary;
			explore_tree(tree, tree_summary);
			explore_tree(rotated_tree, rotated_tree_summary);
			// Saving data
			save_data(tree_summary, output_path + "regular_tree_"+ std::to_string(jittered_iteration) +"/tree_"+iteration_string);
			save_data(rotated_tree_summary, output_path + "rotated_tree_"+std::to_string(jittered_iteration)+"/tree_"+iteration_string);
		}
		/*for(int i = 0; i< simulation_number_of_points; i++){
			Point pt;
			pt = dataset[i];
			for(int j = 0; j<simulation_dimension; j++){
				cout << pt.coordinates[j] << "\t";
			}
			cout << endl;
		}*/
	}
	return 0;	
	
}
