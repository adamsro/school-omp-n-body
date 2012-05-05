import os
#import system

def main():
    for i in [1,2,4,6,8,10,12,14,16]:
        os.system("g++ -Wall -DNUMTHREADS=%d main.cpp -o main -fopenmp -lm" % (i))
        os.system("./main")

if __name__ == "__main__":
    main()
