#include <bits/stdc++.h>
using namespace std;
int main(int argc, char const *argv[])
{
    int length;
    int breadth;
    int grid_x;
    int grid_y;
    cout << "Enter the value of length and breadth of the Lid-Driven Cavity" << endl;
    cin >> length >> breadth;
    cout << "Enter the value for M*N for discretization--(Note :- If you choose Length=Breadth then do M=N)" << endl;
    cin >> grid_x >> grid_y;
    double delta_x;
    double delta_y;
    double Re;
    double U;
    cout << "Enter the value of Reynold's Number" << endl;
    cin >> Re;
    cout << "Enter the value of Free Stream Velocity" << endl;
    cin >> U;
    double mew;
    cout << "Enter Viscocity value (mew)";
    cin >> mew;
    double time_gap;
    cout << "Enter the time gap in the Time-Dependent-Problem" << endl;
    cin >> time_gap;
    delta_x = length / double(grid_x - 1);
    delta_y = breadth / double(grid_y - 1);
    vector<vector<double>> stream(grid_x, vector<double>(grid_y, 0.00));
    vector<vector<double>> vort(grid_x, vector<double>(grid_y, 0.00));

    vector<vector<double>> stream_new(grid_x, vector<double>(grid_y, 0.00));
    vector<vector<double>> vort_new(grid_x, vector<double>(grid_y, 0.00));
    double er;
    cout << "Enter your permissible error for achieving steady-state" << endl;
    cin >> er;
    double er1 = 10;
    double er2 = 10;

    // boundary conditions for the Lid-Driven Cavity:-
    // BC1:-No Slip Walls;-->stream(i,j)=0,vort(i,j)=0 already declared..
    // bC2 :-Moving wall
    for (int j = 0; j < grid_y; j++)
    {
        vort[grid_y - 1][j] = 2 * (stream[grid_y - 1][j] - delta_y - stream[grid_y - 1][j - 1]) / double(pow(delta_y, 2));
    }
    vort_new = vort;

    // Simulation for Interior Points:-
    while (er2 > er || er1 > er)
    {
        
        double den_vort = (1 / double(time_gap)) + (2 / double(Re * pow(delta_x, 2))) + (2 / double(Re * pow(delta_y, 2)));
        for (int i = 1; i < grid_x - 1; i++)
        {
            for (int j = 1; j < grid_y - 1; j++)
            {
                double a = vort[i][j] / double(time_gap);

                double b = (stream_new[i][j+1] - stream_new[i][j-1]) * (vort_new[i+1][j] - vort_new[i-1][j]);
                b = b / double(4 * delta_x * delta_y);

                double c = (stream_new[i+1][j] - stream_new[i-1][j]) * (vort_new[i][j+1] - vort_new[i][j-1]);
                c = c / double(4 * delta_x * delta_y);

                double d = (vort_new[i][j+1] + vort_new[i][j-1]) / double(pow(delta_x, 2));
                d += (vort_new[i+1][j] + vort_new[i-1][j]) / double(pow(delta_y, 2));
                d = d / double(Re);

                vort_new[i][j] = a + b - c + d;
                vort_new[i][j] = vort_new[i][j] / double(den_vort);
            }
        }
        // new vortex fomation at lower wall
        for (int j = 0; j < grid_y; j++)
        {
            vort_new[0][j] = -2 * (stream[1][j] - stream[0][j]) / double(pow(delta_y, 2));
        }

        // new vortex fomation at left wall
        for (int i = 0; i < grid_x; i++)
        {
            vort_new[i][0] = -2 * (stream[i][1] - stream[i][0]) / double(pow(delta_x, 2));
        }

        // new vortex fomation at right wall
        for (int i = 0; i < grid_x; i++)
        {
            vort_new[i][grid_y - 1] = 2 * (stream[i][grid_y - 1] - stream[i][grid_y - 2]) / double(pow(delta_x, 2));
        }
        // new vortex for the upper wall
        for (int j = 0; j < grid_y; j++)
        {
            vort_new[grid_y - 1][j] = 2 * (stream[grid_y - 1][j] - delta_y - stream[grid_y - 2][j]) / double(pow(delta_y, 2));
        }

        double den_stream = (2 / double(pow(delta_x, 2))) + (2 / double(pow(delta_y, 2)));

        for (int i = 1; i < grid_x - 1; i++)
        {
            for (int j = 1; j < grid_y - 1; j++)
            {
                double x = (stream_new[i][j+1] + stream_new[i][j-1]) / double(pow(delta_x, 2));
                x += (stream_new[i+1][j] + stream_new[i-1][j]) / double(pow(delta_y, 2));
                x += vort_new[i][j];
                stream_new[i][j] = x;
                stream_new[i][j] = stream_new[i][j] / double(den_stream);
            }
        }
        er1 = 0, er2 = 0;
        for (int i = 1; i < grid_x-1; i++)
        {
            for (int j = 1; j < grid_y-1; j++)
            {
                if(abs(vort_new[i][j]-vort[i][j])>er1){
                    er1=abs(vort_new[i][j]-vort[i][j]);
                }
                if(abs(stream_new[i][j]-stream[i][j])>er2){
                    er2=abs(stream_new[i][j]-stream[i][j]);
                }
            }
            
        }
        cout <<"Error for Vorticity--->"<< er1 <<"---->  Error for Stream Function--->"<< er2 << endl;
        vort = vort_new;
        stream = stream_new;
    }
    cout << "done1" << endl;

    // calculation for U-direction and V-drection Velocities:-
    vector<vector<double>> U_direct(grid_x, vector<double>(grid_y, 0.00));
    vector<vector<double>> V_direct(grid_x, vector<double>(grid_y, 0.00));
    for (int i = 0; i < grid_x; i++)
    {
        U_direct[grid_x-1][i] = U;
    }

    for (int i = 1; i < grid_x - 1; i++)
    {
        for (int j = 1; j < grid_y - 1; j++)
        {
            U_direct[i][j] = (stream_new[i+1][j] - stream_new[i-1][j]) / double(2 * delta_y);
            U_direct[i][j] = U_direct[i][j] * U;
            V_direct[i][j] = -(stream_new[i][j+1] - stream_new[i][j-1]) / double(2 * delta_x);
            V_direct[i][j] = V_direct[i][j] * U;
        }
    }
    cout << "done2" << endl;

    vector<vector<double>> p(grid_x, vector<double>(grid_y, 0.00));
    double deno = (2 / double(pow(delta_x, 2))) + (2 / double(pow(delta_y, 2)));
    for (int i = 1; i < grid_x - 1; i++)
    {
        for (int j = 1; j < grid_y - 1; j++)
        {
            double psi = pow((stream_new[i + 1][j + 1] - stream_new[i +1][j - 1] - stream_new[i - 1][j + 1] + stream_new[i - 1][j - 1]) / double(4 * delta_x * delta_y), 2);
            double psi2 = ((stream_new[i][j+1] - (2 * stream_new[i][j]) + stream_new[i][j-1]) / double(pow(delta_x, 2))) * ((stream_new[i+1][j] - (2 * stream_new[i][j]) + stream_new[i-1][j]) / double(pow(delta_y, 2)));
            p[i][j] =2 * psi - (2 * psi2);
            p[i][j] = (p[i][j] * mew * Re * U) / length;
        }
    }
    cout << "done3" << endl;

    for (int j = 0; j < grid_y; j++)
    {
        p[0][j] = p[1][j];
        p[grid_y - 1][j] = p[grid_y - 2][j];
    }

    for (int i = 0; i < grid_x; i++)
    {
        p[i][0] = p[i][1];
        p[grid_x - 1][0] = p[grid_x - 2][0];
    }


    // Producting output in file
    //Contour Plots:---------------------
    ofstream out("Assignment_3.plt");
    string var = "VARIABLES=\"X_VAL\",\"Y_VAL\",\"U_VEl\",\"V-VEl\",\"Pressure\",\"Psi\",\"Vorticity\"\n";
    string zone = "ZONE F=POINT\n";
    string col = to_string(grid_x);
    string row = to_string(grid_y);
    string i = "I=" + col + ", J=" + row;
    out << var << zone << i << endl;
    for (int i = 0; i < grid_x; i++)
    {
        for (int j = 0; j < grid_y; j++)
        {
            double xi = i * delta_x;
            double yi = j * delta_y;
            string x = to_string(yi);
            x = x + ' ' + ' ';
            string y = to_string(xi);
            y = y + ' ' + ' ';
            string U_vel = to_string(U_direct[i][j]);
            U_vel = U_vel + ' ' + ' ';
            string V_vel = to_string(V_direct[i][j]);
            V_vel = V_vel + ' ' + ' ';
            string pressure = to_string(p[i][j]);
            pressure = pressure + ' ' + ' ';
            string psi = to_string(stream_new[i][j]);
            psi = psi + ' ' + ' ';
            string vorti = to_string(vort_new[i][j]);
            vorti = vorti + ' ' + ' ';

            out << x << y << U_vel << V_vel << pressure <<psi<<vorti<< endl;
        }
    }
    cout << "done4" << endl;

// U_vs_Y(mid-plane) PLOT-------------------------->
    ofstream outs("U_vs_Y(mid-plane).plt");
    string vars = "VARIABLES=\"Y_VAL\",\"U_VEl\"";
    string zones = "ZONE F=POINT\n";
    string rows = to_string(grid_y);
    string is ="J=" + row;
    outs << vars << zones << is << endl;
    for (int i = 0; i < grid_x; i++)
    {
            double xi = i * delta_x;
            string y = to_string(xi);
            y = y + ' ' + ' ';
            string U_vel = to_string(U_direct[i][50]);
            U_vel = U_vel + ' ' + ' ';
            outs<<y << U_vel << endl;
    }
    cout << "done5" << endl;

    //  v_vs_Y(mid-plane) PLOT-------------------------->
    ofstream outsv("v_vs_Y(mid-plane).plt");
    string varsv = "VARIABLES=\"Y_VAL\",\"V_VEl\"";
    string zonesv = "ZONE F=POINT\n";
    string rowsv = to_string(grid_y);
    string isv ="J=" + row;
    outsv << varsv << zonesv << isv << endl;
    for (int i = 0; i < grid_x; i++)
    {
            double xi = i * delta_x;
            string y = to_string(xi);
            y = y + ' ' + ' '; 
            string U_vel = to_string(V_direct[i][50]);
            U_vel = U_vel + ' ' + ' ';
            outsv<<y << U_vel << endl;
    }
    cout << "done6" << endl;

    return 0;
}
