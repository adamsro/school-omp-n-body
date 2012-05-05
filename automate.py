import os
#import system

def main():
    for i in [1,2,4,8]:
        os.system("g++ -Wall -DNUMTHREADS=%d main.cpp -o main -fopenmp -lm" % (i))
        os.system("./main")

if __name__ == "__main__":
    main()
