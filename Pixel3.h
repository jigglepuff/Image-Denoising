//Karen Wang
//wangkh@usc.edu
//Sept 18 2016

#ifndef Header_h
#define Header_h
#include <vector>
using namespace std;


/* 
This class holds 3 properties of a pixel- 1) x & y position, 
2) RGB intensity and 3) Pixel # or Pixel count (used to compute each pixel's 8 neighbors)
*/

class Pixel{
private:
    int row;
    int col;
    int count;
    int R;
    int G;
    int B;
    
public:
    Pixel(){}

    Pixel(int x,int y)
    {
        row=x;
        col=y;
    }
    
    int getRow()
    {
        return row;
    }
    
    int getCol()
    {
        return col;
    }
    
    int getPixelCount()
    {
        return count;
    }
    
    int getRed()
    {
        return R;
    }
    
    int getGreen()
    {
        return G;
    }
    
    
    int getBlue()
    {
        return B;
    }
    
    void setRow(int inputRow)
    {
        row=inputRow;
    }
    
    void setCol(int inputCol)
    {
        col=inputCol;
    }
    
    void setPixelCount(int PixelCount)
    {
        count=PixelCount;
    }
    
    void setRed(int inputRed)
    {
        R=inputRed;
    }
    
    void setGreen(int inputGreen)
    {
        G=inputGreen;
    }
    
    void setBlue(int inputBlue)
    {
        B=inputBlue;
    }
    
    //Identifying pixel numbers for current pixel's 7 neighbors
    vector<int> getNeighbors(int imageWidth,int windowdim)
    {
        vector<int> neighbors;
        
        for(int row=-windowdim/2; row<windowdim/2+1;row++)
        {
            for(int col=-windowdim/2; col<windowdim/2+1;col++)
            {
                neighbors.push_back(count+(imageWidth*row)+col);
            }
        }
        
        return neighbors;
    }
};

#endif /* Header_h */
