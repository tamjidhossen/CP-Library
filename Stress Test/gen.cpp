#include<bits/stdc++.h>
using namespace std;
int rnd(int a, int b) { // generates random number in range a to b
    return a + rand() % (b - a + 1);
}
int main()
{
    srand(time(NULL)); // commenting this out will result in same random values each time
    cout<<1<<endl; // it's for the number of test cases (should be 1). We will iterate in the bash script.
    // place the input structure of the problem below
    int a = rnd(1, 100), b = (1, 100);
    cout<<a<<" "<<b<<endl; 
    /* 
        Example input:
        1
        20 2
    */
    return 0;
}