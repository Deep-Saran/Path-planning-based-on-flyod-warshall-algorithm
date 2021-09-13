#include "network.h"
#include <cmath>
#include <string>
using namespace std;

#define network_size 303
#define INF 9999999 

#define conditional_swap(a,b,c,d) if(b<a){ a=b;c=d; } //for sawpping matrix elemnts

class Path{
    public:
        Path(string a, string b); // constructor
        string intial_charger;
        string end_charger;
        //const int network_size = 303;
        int full_charge = 320;
        int constant_speed = 105;
        int start = 0;
        int end = 0;
        int u_limit = INF;
        int l_limit = -1;
        double radius = 6356.752; //in Km
        double a[network_size][network_size], A[network_size][network_size],
                dist_2D[network_size][network_size], ratio_S2R[network_size];         

        bool isValid();
        double degree_to_radian(double degree);
        double distance(double lat1, double lat2, double lon1, double lon2);
        double speed_to_rate(double charge_rate);
        // void conditional_swap(int &a, int &b, int &c, int &d);
        string output_format(int i, int j);
        void flyod_Warshall();

};

Path::Path(string a, string b){// Defining the constructor 
    intial_charger = a;
    end_charger = b;
    isValid();
}

bool Path::isValid(){
    //matching given input to the network
    while (start < network.size() && network[start].name != intial_charger){ ++start;}
    while (end < network.size() && network[end].name != end_charger){ ++end;}
    if(start <= network.size() || end>= network.size()){
        return false;
    }
    return true;
    // cout << start << endl;
    // cout << end << endl;
}

double Path::degree_to_radian(double degree){ //converting lattitudes and longitudes from degrees to radians
    return (degree/180)*3.1415926535;
}

double Path::distance(double lat1, double lat2, double lon1, double lon2){ // distance between two points on earth
    double ans = radius*(2*asin(sqrt(pow(((degree_to_radian(lat1)-degree_to_radian(lat2))/2),2)
                + cos(degree_to_radian(lat1))*cos(degree_to_radian(lat2)) * 
                pow(((degree_to_radian(lon1)-degree_to_radian(lon2))/2),2))));
    return ans;
}

double Path::speed_to_rate(double charge_rate){
    return constant_speed/charge_rate;
}

// void Path::conditional_swap(int &a, int &b, int &c, int &d){ if(b<a){ a=b;c=d; }}

//flyod_warshall algorithm for optimal path (ALL pair shortest path)
void Path::flyod_Warshall(){ // input distance array and ratio
    for(int i =0;i < network_size;++i){
        ratio_S2R[i] = speed_to_rate(network[i].rate);
        for(int j=0; j< network_size; ++j){ //Distance between every pair of node
            dist_2D[i][j] = distance(network[i].lat,network[j].lat,network[i].lon,network[j].lon);
            if(dist_2D[i][j] > full_charge) { dist_2D[i][j] = INF;}
            // cout << dist_2D[i][j] << endl;
        }
        // cout << ratio_S2R[i] << endl;
    }
    //INF for absense of edge and -1 for delf loop
    for(int i=0; i< network_size; ++i){
        for(int j=0; j< network_size; ++j){
            a[i][j] = INF;
            A[i][j] = -1;
            a[i][j] = dist_2D[i][j] * (ratio_S2R[i]+1.); 
            //cout << a[i][j] << endl;
            //cout << A[i][j] << endl;
            }
        }

    /*Multiple cases for finding the shortest time path between the nodes: 
    -we reach destination with zero range and never charge to full range
    -we start with full range, reach with zero range, and never reach zero range or full range in between nodes 
     (computing all possibilities based on path length between nodes)
    -we start with full range, reach with zero range, and never reach full range in between
    -we start with full range and reach destination with full range.*/
    for(int k=0; k<network_size;++k){
        for(int i =0;i < network_size;++i){
            for(int j=0; j< network_size; ++j){
                conditional_swap(a[i][j], a[i][k] + a[k][j], A[i][j], k); 
                if (dist_2D[i][k] + dist_2D[k][j] >= full_charge){
                    conditional_swap(a[i][j], full_charge +  
                                    (dist_2D[i][k]+dist_2D[k][j]-full_charge)*(ratio_S2R[k]+1.), A[i][j], k);
                }
                conditional_swap(a[i][j], a[i][k] + a[k][j], A[i][j], k);
                a[i][j] = a[i][j] + ratio_S2R[j] * full_charge;
                //cout << a[i][j] << endl;
                conditional_swap(a[i][j], dist_2D[i][j] * (1. + ratio_S2R[j]), A[i][j], -2);
                conditional_swap(a[i][j], a[i][k] + a[k][j], A[i][j], k);
            }
        }
    }

    for(int i =0;i < network_size;++i){
        conditional_swap(u_limit, a[start][i] + a[i][end], l_limit, i);
    }
    // cout << l_limit << endl;
    cout << intial_charger + output_format(start, l_limit) + output_format(l_limit, end) + 
            ", " + end_charger << endl;
}

string Path::output_format(int i, int j) {
    if (A[i][j] == -1)
        return output_format(i, j) + ", " + to_string(ratio_S2R[j]/constant_speed);
    if (A[i][j] == -2)
        return ", " + network[j].name + ", " + to_string(dist_2D[i][j] * ratio_S2R[j] / constant_speed);
    return output_format(i, A[i][j]) + output_format(A[i][j], j);
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "Error: requires initial and final supercharger names" << endl;        
        return -1;
    }
    
    string initial_charger_name = argv[1];
    string goal_charger_name = argv[2];

    Path p(initial_charger_name,goal_charger_name);
    p.flyod_Warshall();

    return 0;
}
