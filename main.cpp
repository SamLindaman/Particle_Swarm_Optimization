/*
 Particle Swarm Optimization
 Beijing Jiaotong DaXue C++ Final Project
 Samuel David Lindaman
 18309004

 Objective:
 Solve several optimization problems
    Min 5 problems
    At most 15 problems
    Solve more problems, get higher score
 Required materials
    Document - analysis of your codes, solved problems, Figures of final results, Comments about PSO
    PPT files
    Source code
 Deadline
    18:00 on 11/6
 
 
 Equations are listed F1-F15 in order of complexity
 PSO code is after the equations
*/

#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#define rand_01 ((float)rand() / (float)RAND_MAX)

// set dimension and number of particles
const int numOfDims = 10;
const int numOfParticles = 20;

using namespace std;



/* ______________________________ Function 1 ______________________________ */
// Sphere function Equation
void F1(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles); // create space fitnesses[i] array
    for(int i = 0; i < numOfParticles; i++){
        for(int j = 0; j < numOfDims; j++){
            t1 = X[i][j];
            t1 *= t1;
            fitnesses[i] += t1;
        }
    }
    p= fitnesses[numOfParticles-1];
}

/* ______________________________ Function 2 ______________________________ */
// Schwefel's 2.22 function Equation
void F2(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1,t2 = 0.0;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                t1 = abs(X[i][j]);
                t2 *= abs(X[i][j]);
                fitnesses[i] += t1 ;//+ t2;
        }
    }
        p= fitnesses[numOfParticles-1];
}

/* ______________________________ Function 3 ______________________________ */

// Schwefel's 1.20 function Equation
void F3(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                t1 = X[i][j-1];
                t1 += t1;
                t1 *= t1;
                fitnesses[i] += t1;
        }
    }
    p= fitnesses[numOfParticles-1];
}

/* ______________________________ Function 4 ______________________________ */

// Schwefel's 2.21 function Equation
void F4(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                t1 = abs(X[i][j]);
                t1 *= t1;
                fitnesses[i] = t1;
        }
    }
    p= fitnesses[numOfParticles-1];
}

/* ______________________________ Function 5 ______________________________ */
// Rosenbrock function Equation
void F5(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float x1, x2, t1, t2;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
    for(int i = 0; i < numOfParticles; i++)
        for(int j = 0; j < numOfDims - 1; j++){
            x1 = X[i][j];
            x2 = X[i][j+1];
            t1 = (x2 - x1 * x1);
            t1 *= t1;
            t1 *= 100;
            t2 = x1 - 1;
            t2 *= t2;
            fitnesses[i] = t1 + t2;
        }
        p= fitnesses[numOfParticles-1];
}

/* ______________________________ Function 6 ______________________________ */
// Step function Equation
void F6(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                t1 = X[i][j] + 0.5;
                t1 *= t1;
                fitnesses[i] += t1;
        }
    }
    p= fitnesses[numOfParticles-1];
}


/* ______________________________ Function 7 ______________________________ */
// Quartic noise function Equation
void F7(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float x1,t1 = 0.0,t2;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                x1 = X[i][j];
                x1 *= x1;
                x1 *= x1;
                x1 *= j;
                t1 += x1;
                t2 = rand() % 2;
                fitnesses[i] = t1 + t2;
                
        }
    }
}

/* ______________________________ Function 8 ______________________________ */
// Schwefel function Equation
void F8(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                fitnesses[i] = ((- X[i][j])* sin (sqrt(abs( X[i][j] ))));
        }
    }
}

/* ______________________________ Function 9 ______________________________ */
// Rastrigin function Equation
void F9(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    double pi = acos(-1);
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                fitnesses[i] += (pow(X[i][j],2) - 10 * cos (2 * pi * X[i][j]) + 10);
        }
    }
}


/* ______________________________ Function 10 ______________________________ */
// Ackley function Equation
void F10(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float x1,t1,t2;
    double pi = acos(-1);
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                x1 = X[i][j];
                x1 *= x1;
                t1 = x1 / j;
                t2 = cos(2 * pi * X[i][j]) / j;
                fitnesses[i] = -20 * exp(-0.2 * sqrt(t1)) - exp(t2) + 20 + exp(1);
        }
    }
}

/* ______________________________ Function 11 ______________________________ */
// Griewank function Equation
void F11(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                t1 = 1;
                fitnesses[i] = pow(X[i][j],2) / 4000 + 1;
                for(int k = 0; k < numOfDims; k++){
                    t1 = t1 * cos(X[i][k] / sqrt(k+1));
                }
                fitnesses[i] = fitnesses[i] - t1;
        }
    }
}


/* ______________________________ Function 12 ______________________________ */
void F12(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float y1,y2,y3,y4,t1,t2,t3,t4;
    float z = 0.0;
    double pi = acos(-1);
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims-1; j++){
                y1 = X[i][1];
                t1 = pow(sin(pi * y1),2);
                t1 *= 10;
                y2 = X[i][j];
                t2 = pow((y2 - 1),2);
                y3 = X[i][j+1];
                t3 = pow(sin ( pi * y3),2);
                t3 *= 10;
                y4 = X[i][numOfDims-1];
                t4 = pow((y4 - 1),2);
                fitnesses[i] = t2 * ( 1 + t3 + t4);
        
                for (int h = 0; h < numOfDims; h++){
                    if(X[i][h] > 10){
                        z = 100 * (X[i][h]-10) * (X[i][h]-10) * (X[i][h]-10) * (X[i][h]-10);
                    } else if(-10 <= X[i][h] && X[i][h] <= 10){
                        z = 0;
                    } else {
                        z = 100 * (- X[i][h] - 10) * (- X[i][h] - 10) * (- X[i][h] - 10) * (- X[i][h] - 10);
                    }
                }
                fitnesses[i] = ((pi / numOfDims)* (t1 * fitnesses[i]))+ z;
            }
                
    }
}


/* ______________________________ Function 13 ______________________________ */
void F13(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float x1,x2,t1,t2,t3,t4,t5,t6;
    float z = 0.0;
    double pi = acos(-1);
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
                x1 = X[i][1];
                t1 = sin ( 3 * pi * x1);
                t1 *= t1;
            for(int j = 0; j < numOfDims; j++){
                x2 = X[i][j];
                t2 = x2 - 1;
                t2 *= t2;
                t3 = sin ( 3 * pi * x2 + 1);
                t3 *= t3;
                
                fitnesses[i] += t2 * (1 + t3);
                
                t4 = x2 - 1;
                t4 *= t4;
                t5 = sin ( 2 * pi * x2);
                t5 *= t5;
                t6 = t4 * t5;
                
                for (int h = 0; h < numOfDims; h++){
                    if(X[i][h] > 5){
                        z = 100 * (X[i][h]-5) * (X[i][h]-5) * (X[i][h]-5) * (X[i][h]-5);
                    } else if(-5 <= X[i][h] && X[i][h] <= 5){
                        z = 0;
                    } else {
                        z = 100 * (- X[i][h] - 5) * (- X[i][h] - 5) * (- X[i][h] - 5) * (- X[i][h] - 5);
                    }
                }
        fitnesses[i] = 0.1 * (t1 + fitnesses[i] + t6) + z;
        }

    }
}


/* ______________________________ Function 14 ______________________________ */
void F14(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float x1,a1,t1;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < 25; j++){
                for(int k = 0; k < 2; k++){
                    x1 = X[i][k];
                    a1 = X[k][j];
                    t1 = x1 - a1;
                    t1 = t1 * t1 * t1 * t1 * t1 * t1;
                    fitnesses[i] += t1;
                    fitnesses[i] = j + fitnesses[i];
                    fitnesses[i] = 1 / fitnesses[i];
                    }
            fitnesses[i] += fitnesses[i];
        }
        fitnesses[i] = 1/500 + fitnesses[i];
        fitnesses[i] = 1 / fitnesses[i];
    }
}


/* ______________________________ Function 15 ______________________________ */
void F15(float X[numOfParticles][numOfDims], float fitnesses[numOfParticles], float &p){
    float a1,b1,t1,t2,t3,t4;
    memset(fitnesses, 0, sizeof (float) * numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < 11; j++){
                a1 = X[1][j];
                b1 = X[2][j];
                t1 = X[i][j] * (b1 * b1 + b1 * X[i][2]);
                t2 = (b1 * b1 + b1 * X[i][3] + X[i][4]);
                t3 = t1 / t2;
                t4 = a1 - t3;
                t4 *= t4;
                fitnesses[i] = t4;
        }
    }

}


/* ____________________________ Particle Swarm Optimization ____________________________ */
// PSO explanation in PPT slides in github repo

void PSO(int numOfIterations, float c1, float c2,
              float xMin[numOfDims], float xMax[numOfDims], float initialPop[numOfParticles][numOfDims],
                float bestArray[], float *gBestFit, float gBest[numOfDims], float* p, int select){
    
    // STEP 1. value declaration
    // V = velocity
    // p = position
    // x = particle coordinates
    float V[numOfParticles][numOfDims] = {0};
    float X[numOfParticles][numOfDims];
    float Vmax[numOfDims];
    float Vmin[numOfDims];
    float pbests[numOfParticles][numOfDims];
    float pbestfits[numOfParticles];
    float fitnesses[numOfParticles];
    float w;
    float minfit;
    int   minfitidx;
    
    //copies initial random values array to initial positions of the particles
    // SETS INITIAL PARTICLE POSITIONS FOR SELECTED EQUATION
    memcpy(X, initialPop, sizeof(float) * numOfParticles * numOfDims);
    
        switch(select){
        case 1:
            F1(X, fitnesses, *p);
            break;
        case 2:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 10;
                xMin[i] = -10;
                }
            F2(X, fitnesses, *p);
            break;
        case 3:
            F3(X, fitnesses, *p);
            break;
        case 4:
            F4(X, fitnesses, *p);
            break;
        case 5:
            F5(X, fitnesses, *p);
            break;
        case 6:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 1.28;
                xMin[i] = -1.28;
                }
            F6(X, fitnesses, *p);
            break;
        case 7:
            F7(X, fitnesses, *p);
            break;
        case 8:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 500;
                xMin[i] = -500;
                }
            F8(X, fitnesses, *p);
            break;
        case 9:
            F9(X, fitnesses, *p);
            break;
        case 10:
            F10(X, fitnesses, *p);
            break;
        case 11:
            F11(X, fitnesses, *p);
            break;
        case 12:
            F12(X, fitnesses, *p);
            break;
        case 13:
            F13(X, fitnesses, *p);
            break;
        case 14:
            F14(X, fitnesses, *p);
            break;
        case 15:
            F15(X, fitnesses, *p);
            break;
        default:
            ;
        }
    
    
    minfit = *min_element(fitnesses, fitnesses + numOfParticles);
    minfitidx = min_element(fitnesses, fitnesses + numOfParticles) - fitnesses;
    
    *gBestFit = minfit; // initialize the global 'best fit' for the particle
    memcpy(gBest, X[minfitidx], sizeof(float) * numOfDims);
    
    // STEP 2. Set velocity limit
    for(int i = 0; i < numOfDims; i++){
        Vmax[i] = 0.2 * (xMax[i] - xMin[i]);
        Vmin[i] = -Vmax[i];
    }

    // STEP 3. Calculate the minimum value of individual history
    // loop iterates 1000 times through the pso algorithm
    for(int t = 0; t < 1000; t++){
        w = 0.2 * t / numOfIterations;

        // STEP 4. set pbest - The best position of the particle
        for(int i = 0; i < numOfParticles; i++){
            if(fitnesses[i] < pbestfits[i]){
                pbestfits[i] = fitnesses[i];
                memcpy(pbests[i], X[i], sizeof(float) * numOfDims);
            }
        }
       
        // STEP 5. Particle update rule
        for(int i = 0; i < numOfParticles; i++){
            for(int j = 0; j < numOfDims; j++){
                
                // velocity
                V[i][j] = min(max((w * V[i][j] + rand_01 * c1 * (pbests[i][j] - X[i][j]) + rand_01 * c2 * (gBest[j] - X[i][j])), Vmin[j]), Vmax[j]);
                
                // position : p = p + v
                X[i][j] = min(max((X[i][j] + V[i][j]), xMin[j]), xMax[j]);
            }
        }

        // GET FINAL EQUATION RESULTS AFTER PSO
        switch(select){
        case 1:
            F1(X, fitnesses, *p);
            break;
        case 2:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 10;
                xMin[i] = -10;
                }
            F2(X, fitnesses, *p);
            break;
        case 3:
            F3(X, fitnesses, *p);
            break;
        case 4:
            F4(X, fitnesses, *p);
            break;
        case 5:
            F5(X, fitnesses, *p);
            break;
        case 6:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 1.28;
                xMin[i] = -1.28;
                }
            F6(X, fitnesses, *p);
            break;
        case 7:
            F7(X, fitnesses, *p);
            break;
        case 8:
            for (int i = 0; i < numOfParticles; i++)
                {
                xMax[i] = 500;
                xMin[i] = -500;
                }
            F8(X, fitnesses, *p);
            break;
        case 9:
            F9(X, fitnesses, *p);
            break;
        case 10:
            F10(X, fitnesses, *p);
            break;
        case 11:
            F11(X, fitnesses, *p);
            break;
        case 12:
            F12(X, fitnesses, *p);
            break;
        case 13:
            F13(X, fitnesses, *p);
            break;
        case 14:
            F14(X, fitnesses, *p);
            break;
        case 15:
            F15(X, fitnesses, *p);
            break;
        default:
            ;
        }

    minfit = *min_element(fitnesses, fitnesses + numOfParticles);
    minfitidx = min_element(fitnesses, fitnesses + numOfParticles) - fitnesses;
    
        // STEP 6. set gbest : The best position of the swarm
        if(minfit < *gBestFit){
            *gBestFit = minfit;
            memcpy(gBest, X[minfitidx], sizeof(float) * numOfDims);
            cout << t << " times fitness is " << minfit << endl;
        }
        
        bestArray[t] = *gBestFit;

    }
}


/* ______________________________ MAIN ______________________________ */

int main(){
    
    time_t t;
    srand((unsigned) time(&t));

    float xMin[numOfParticles], xMax[numOfParticles]; // value declaration
    float initPop[numOfParticles][numOfDims];
    float bests[1000];
    float gBestFit, point;
    float gBest[numOfDims];
    int select;
    
    // set lower and upper value
    for(int i = 0; i < 10; i++){
        xMax[i] = 100;
        xMin[i] = -100;
        }
        
    //set all particles to random values
    for(int i = 0; i < 20; i++)
        for(int j = 0; j < 10; j++){
               initPop[i][j] = rand() % (100 + 100 + 1) - 100;
        }
    
    
    while(1){ // Select the Questions.
      cout << endl << "\n\n\n\n\nPSO         " << endl;
      cout << endl << "=======================================================" << endl << endl;
      cout << "01. Question 1 - Sphere function Equation" << endl;
      cout << "02. Question 2 - Schwefel's 2.22 function Equation" << endl;
      cout << "03. Question 3 - Schwefel's 1.20 function Equation" << endl;
      cout << "04. Question 4 - Schwefel's 2.21 function Equation" << endl;
      cout << "05. Question 5 - Rosenbrock function Equation" << endl;
      cout << "06. Question 6 - Step function Equation" << endl;
      cout << "07. Question 7 - Quartic noise function Equation" << endl;
      cout << "08. Question 8 - Schwefel function Equation" << endl;
      cout << "09. Question 9 - Rastrigin function Equation" << endl;
      cout << "10. Question 10 - Ackley function Equation" << endl;
      cout << "11. Question 11 - Griewank function Equation" << endl;
      cout << "12. Question 12" << endl;
      cout << "13. Question 13" << endl;
      cout << "14. Question 14" << endl;
      cout << "15. Question 15" << endl;
      cout << "Enter any other number to exit" << endl << endl;
      cout << "=======================================================" << endl;
      cout << "Command : ";
      
      cin >> select;
      cout << endl;
          if(select>15 || select<=0){
            cout<<"Goodbye!\n";
            exit(0);
        }else{
            PSO(1000, 2, 2, xMin, xMax, initPop, bests, &gBestFit, gBest, &point, select);
            cout << endl << "Result : " << gBestFit << endl;
              cout << endl << "================== Question." << select << " Finish ==================" << endl;
        }
    }
    return 0;
}
