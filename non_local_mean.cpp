//Karen Wang
//wangkh@usc.edu
//Sept 18 2016

#include "Pixel3.h"
#include <iostream>
#include <cmath>



/*
Function: computeNeighborhoodMean()
Input: Center pixel (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), local neighborhood dimension (int)
Purpose: Calculates mean intensity (for R, G or B specified by input) of a local neighborhood specified by input neighborhoodDim
Output: integer that represents mean intensity of local neighborhood
*/

int computeNeighborhoodMean(Pixel center, Pixel* PixArray, int color, int height, int width, int neighborhoodDim)
{
    int average=0;
    int row =center.getRow();
    int col = center.getCol();
    int neighborhoodSize=neighborhoodDim*neighborhoodDim;
    int bounds=neighborhoodDim/2;
    vector<int> neighbors;
    
    //cout<<"center=("<<center.getRow()<<","<<center.getCol()<<")"<<endl;
    
    if (row<bounds || col<bounds || row>height-1-bounds || col>width-1-bounds) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                average=(center.getRed())*neighborhoodSize;
                break;
            case 1:
                average=(center.getGreen())*neighborhoodSize;
                break;
            case 2:
                average=(center.getBlue())*neighborhoodSize;
                break;
                
            default:
                break;
        }
    }
    else//within bounds
    {
        neighbors=center.getNeighbors(width,neighborhoodDim);
        for (int i=0; i<neighborhoodSize; i++) {
            switch (color) {
                case 0:
                    average+=PixArray[neighbors[i]].getRed();
                    if (row==10 && col==10) {
                        //cout<<PixArray[neighbors[i]].getRed()<<endl;
                    }
                    break;
                case 1:
                    average+=PixArray[neighbors[i]].getGreen();
                    break;
                    
                case 2:
                    average+=PixArray[neighbors[i]].getBlue();
                    break;
                    
                default:
                    break;
            }
        }//end of neighbours loop
        for (int i=0; i<9; i++) {
            neighbors.pop_back();
        }
    }
    
    //Uncomment to debug: 
    //cout<<"mean= "<<average/neighborhoodSize<<endl;
    return average/neighborhoodSize;
}


/*
Function: nonLocalMean()
Input: Center pixel (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), non-local window dimension (int),
       Gaussian filter variance (float), local neighborhood dimension (int)
Purpose: Calls computeNeighborhoodMean(), and computes weighted average of all the local neighborhoods specified by user input
Output: integer that represents the weighted average of the local neighborhood means.
*/
int NonLocalMean(Pixel center, Pixel* PixArray, int color, int height, int width, int windowDim, float variance, int neighborhoodDim)
{
    float weightedAvg=0;
    int row =center.getRow();
    int col = center.getCol();
    int bounds=windowDim/2;
    int windowSize=windowDim*windowDim;
    
    float gaussianRatio=0;
    float exponent=0;
    float gaussianSum=0;
    int currLocalMean=0, neighborLocalMean=0;
    
    vector<int> neighbors;
    vector<float> weights;

    //Uncomment to debug: 
    /*
    cout<<"current=("<<row<<","<<col<<","<<color<<")"<<endl;
    cout<<"\n\ni=("<<row<<","<<col<<")"<<endl;
    */
    if (row<bounds || col<bounds || row>height-1-bounds || col>width-1-bounds) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                //Uncomment to debug: 
                //weightedAvg=(center.getRed());
                weightedAvg=computeNeighborhoodMean(PixArray[center.getPixelCount()],PixArray,color,height,width,neighborhoodDim);
                break;
            case 1:
                 //Uncomment to debug: 
                //weightedAvg=(center.getGreen());
                weightedAvg=computeNeighborhoodMean(PixArray[center.getPixelCount()],PixArray,color,height,width,neighborhoodDim);
                break;
            case 2:
                //Uncomment to debug: 
                //weightedAvg=(center.getBlue());
                weightedAvg=computeNeighborhoodMean(PixArray[center.getPixelCount()],PixArray,color,height,width,neighborhoodDim);
                break;
                
            default:
                break;
        }
    }
    else //within bounds
    {
        //Computes center's local mean
        currLocalMean=computeNeighborhoodMean(PixArray[center.getPixelCount()],PixArray,color,height,width,neighborhoodDim);
        
        //Uncomment to debug: 
        /*
        if (center.getRow()==10 && center.getCol()==10) {
            cout<<"current=("<<center.getRow()<<","<<center.getCol()<<")"<<endl;
            cout<<"currLocalMean="<<currLocalMean<<endl;
        }
        */

        
        //Stores all neighbors' pixel #;
        neighbors=center.getNeighbors(width,windowDim);
        
        //Computes all neighbors' local mean
        gaussianSum=0;
        exponent=0;
        gaussianRatio=0;
        for (int i=0; i<windowSize; i++) {
            neighborLocalMean=computeNeighborhoodMean(PixArray[neighbors.at(i)],PixArray,color,height,width,neighborhoodDim);
            exponent=-abs(neighborLocalMean-currLocalMean)/(pow(10*variance, 2));
            gaussianRatio=(1/(2*3.141593*pow(variance, 2)))*exp(exponent);
            weights.push_back(gaussianRatio);
            
            //Uncomment to debug: 
            /*
            if (center.getRow()==10 && center.getCol()==10) {
                cout<<"neighborLocalMean="<<neighborLocalMean<<endl;
                cout<<neighborLocalMean<<endl;
                cout<<"gaussianRatio @"<<i<<"="<<gaussianRatio<<endl;
                cout<<gaussianRatio<<endl;
            }
            */
            gaussianSum+=gaussianRatio;
        }//end of for loop

        //Uncomment to debug: 
        /*
        if (center.getRow()==10 && center.getCol()==10) {
            cout<<"gaussianSum="<<gaussianSum<<endl;
        }
        */
    
        //Computes weighted average from Gaussian weights;
        weightedAvg=0.0;
        for (int i=0; i<windowSize; i++) {
            switch (color) {
                case 0:
                    weightedAvg+=(weights.at(i)/gaussianSum)*float(PixArray[neighbors.at(i)].getRed());
                    //Uncomment to debug: 
                    /*
                    if (center.getRow()==10 && center.getCol()==10) {
                        cout<<"weighted gaussian@"<<i<<"="<<(weights.at(i)/gaussianSum)<<endl;
                        cout<<"unfiltered intensity="<<float(PixArray[neighbors.at(i)].getRed())<<endl;
                        cout<<"multiply="<<(weights.at(i)/gaussianSum)*float(PixArray[neighbors.at(i)].getRed())<<endl;
                        cout<<"weightedAvg="<<weightedAvg<<endl;
                    }
                    */
                break;
                    
                case 1:
                    weightedAvg+=(weights.at(i)/gaussianSum)*float(PixArray[neighbors.at(i)].getGreen());
                    break;
                    
                case 2:
                    weightedAvg+=(weights.at(i)/gaussianSum)*float(PixArray[neighbors.at(i)].getBlue());
                    break;
                default:
                    break;
            }//end of switch
        }//end of for loop
    }//end of else case
    
    weightedAvg=int(weightedAvg);
    //Uncomment to debug: 
    /*
    if (center.getRow()==10 && center.getCol()==10) {
        cout<<"weightedAvg="<<weightedAvg<<endl;
        cout<<weightedAvg<<endl;
    }
    */
    
    for (int i=0; i<windowSize; i++) {
        neighbors.pop_back();
        weights.pop_back();
    }
    

    return weightedAvg;
}



/////////////////////// main() ///////////////////////

int main(int argc, const char * argv[]) {
    
    // Define file pointer and variables
    FILE *file;
    int BytesPerPixel;//greyscale or RGB
    int SizeWidth = 256; //default width
    int SizeHeight= 256; //default height
    
    // Check for proper syntax
    if (argc < 3){
        cout << "Syntax Error - Incorrect Parameter Usage:" << endl;
        cout << "program_name [input_image.raw] [output_image.raw] [BytesPerPixel] [Width] [Height]" << endl;
        return 0;
    }
    
    // Check if image is grayscale or color
    if (argc < 4){
        BytesPerPixel = 1; // default is grey image
    }
    else {
        BytesPerPixel = atoi(argv[3]);
        // Check if size is specified
        if (argc >= 5){
            SizeWidth = atoi(argv[4]);
            SizeHeight= atoi(argv[5]);
        }
    }
    
    // Allocate original image data array
    unsigned char Imagedata[SizeHeight][SizeWidth][BytesPerPixel];
    unsigned char originaldata[SizeHeight][SizeWidth][BytesPerPixel];
    
    // Read image (filename specified by first argument) into image data matrix
    if (!(file=fopen(argv[1],"rb"))) {
        cout << "Cannot open file: " << argv[1] <<endl;
        exit(1);
    }
    fread(Imagedata, sizeof(unsigned char), SizeHeight*SizeWidth*BytesPerPixel, file);
    fclose(file);
    
    
    // Read image (filename specified by first argument) into image data matrix
    if (!(file=fopen(argv[6],"rb"))) {
        cout << "Cannot open file: " << argv[1] <<endl;
        exit(1);
    }
    fread(originaldata, sizeof(unsigned char), SizeHeight*SizeWidth*BytesPerPixel, file);
    fclose(file);
    


    //Variable declariation
    int newWidth;
    int newHeight;
    
    if (argc<6) {
        cout << "Syntax Error - Incorrect Parameter Usage:" << endl;
        cout << "program_name input_image.raw output_image.raw [BytesPerPixel] [InputImageWidth] [InputImageHeight]" << endl;
        exit(1);
    }
    
    //User inputs:
    cout<<"User Inputs:"<<endl;
    //cout<<argv[0]<<endl;
    cout<<"Input filename: "<<argv[1]<<endl;
    cout<<"Output filename: "<<argv[2]<<endl;
    cout<<"Bytes Per Pixel: "<<argv[3]<<endl;
    cout<<"Input Image Width: "<<argv[4]<<endl;
    cout<<"Input Image Height: "<<argv[5]<<endl;
    newHeight=SizeHeight;
    newWidth=SizeWidth;
    
    //Processing starts here:
    //Initializing Output Image to all black pixels
    int pixelCount=0;
    Pixel* PixelArray;
    PixelArray = new Pixel [SizeWidth*SizeHeight];
    //Pixel PixelArray[262144];
    unsigned char tempImage[newHeight][newWidth][BytesPerPixel];
    for (int i=0; i<newHeight; i++)
    {
        for(int j=0; j<newWidth; j++)
        {
            PixelArray[pixelCount].setRow(i);
            PixelArray[pixelCount].setCol(j);
            PixelArray[pixelCount].setPixelCount(pixelCount);
            PixelArray[pixelCount].setRed(int(Imagedata[i][j][0]));
            PixelArray[pixelCount].setGreen(int(Imagedata[i][j][1]));
            PixelArray[pixelCount].setBlue(int(Imagedata[i][j][2]));
            for (int k=0; k<BytesPerPixel; k++) {
                //tempImage[i][j][k]=Imagedata[i][j][k];
                tempImage[i][j][k]=0;
            }
            pixelCount++;
        }
    }
    
    //Move window through image
    pixelCount=0;
    float PSNR_R=0,PSNR_G=0,PSNR_B=0;
    float MSE_R=0,MSE_G=0,MSE_B=0;
    //SSIM variables
    float muX_R=0, muY_R=0,varX_R=0, varY_R=0, varXY_R=0;
    float muX_G=0, muY_G=0,varX_G=0, varY_G=0, varXY_G=0;
    float muX_B=0, muY_B=0,varX_B=0, varY_B=0, varXY_B=0;
    
    for (int i=0; i<newHeight; i++)
    {
        for(int j=0; j<newWidth; j++)
        {
            for (int k=0; k<BytesPerPixel; k++) {
                //Linear 3x3
                //tempImage[i][j][k]=char(LinearFilter(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth));
                
                //Local Mean 7x7
                //tempImage[i][j][k]=computeNeighborhoodMean(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,7);

                //Local Mean 5x5
                tempImage[i][j][k]=NonLocalMean(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth, 15, 0.25, 5);

            }
            
            //MSE Calculation
            MSE_R+=pow(tempImage[i][j][0]-originaldata[i][j][0],2);
            MSE_G+=pow(tempImage[i][j][1]-originaldata[i][j][1],2);
            MSE_B+=pow(tempImage[i][j][2]-originaldata[i][j][2],2);
            
            //SSIM Index
            muX_R+=originaldata[i][j][0];
            muX_G+=originaldata[i][j][1];
            muX_B+=originaldata[i][j][2];
            
            muY_R+=tempImage[i][j][0];
            muY_G+=tempImage[i][j][1];
            muY_B+=tempImage[i][j][2];
            
            pixelCount++;
        }
    }
    
    MSE_R=MSE_R/float(SizeHeight*SizeWidth);
    MSE_G=MSE_G/float(SizeHeight*SizeWidth);
    MSE_B=MSE_B/float(SizeHeight*SizeWidth);
    
    PSNR_R=10*log10((255.0*255.0)/MSE_R);
    PSNR_G=10*log10((255.0*255.0)/MSE_G);
    PSNR_B=10*log10((255.0*255.0)/MSE_B);

    //SSIM muX, muY calculation
    muX_R=muX_R/float(SizeHeight*SizeWidth);//original image av
    muX_G=muX_G/float(SizeHeight*SizeWidth);
    muX_B=muX_B/float(SizeHeight*SizeWidth);
    
    muY_R=muY_R/float(SizeHeight*SizeWidth);//distorted image av
    muY_G=muY_G/float(SizeHeight*SizeWidth);
    muY_B=muY_B/float(SizeHeight*SizeWidth);
    
    //SSIM Index
    float SSIM_R=0, SSIM_G=0, SSIM_B=0;
    float L=pow(2, 8)-1; //2^(#ofbits per pixel)-1
    float k1=0.01;
    float k2=0.03;
    float c1=pow(k1*L,2);
    float c2=pow(k2*L,2);
    
    for (int i=0; i<newHeight; i++)
    {
        for(int j=0; j<newWidth; j++)
        {
            for (int k=0; k<BytesPerPixel; k++)
            {
                varX_R+=pow(float(originaldata[i][j][0])-muX_R, 2.0);// variance of orignal image
                varX_G+=pow(float(originaldata[i][j][1])-muX_G, 2.0);
                varX_B+=pow(float(originaldata[i][j][2])-muX_B, 2.0);
                
                varY_R+=pow(float(tempImage[i][j][0])-muY_R, 2.0);// variance of distorted image
                varY_G+=pow(float(tempImage[i][j][1])-muY_G, 2.0);
                varY_B+=pow(float(tempImage[i][j][2])-muY_B, 2.0);
                
                varXY_R+=(float(originaldata[i][j][0]-muX_R)*float(tempImage[i][j][0]-muY_R));
                varXY_G+=(float(originaldata[i][j][1]-muX_G)*float(tempImage[i][j][1]-muY_G));
                varXY_B+=(float(originaldata[i][j][2]-muX_B)*float(tempImage[i][j][2]-muY_B));
            }
        }
    }
    
    varX_R=varX_R/float(SizeHeight*SizeWidth);
    varX_G=varX_G/float(SizeHeight*SizeWidth);
    varX_B=varX_B/float(SizeHeight*SizeWidth);
    
    varY_R=varY_R/float(SizeHeight*SizeWidth);
    varY_G=varY_G/float(SizeHeight*SizeWidth);
    varY_B=varY_B/float(SizeHeight*SizeWidth);
    
    varXY_R=varXY_R/float(pow(SizeHeight,2)*pow(SizeWidth,2));
    varXY_G=varXY_G/float(pow(SizeHeight,2)*pow(SizeWidth,2));
    varXY_B=varXY_B/float(pow(SizeHeight,2)*pow(SizeWidth,2));
    
    SSIM_R=((2*muX_R*muY_R+c1)*(2*varXY_R+c2))/((pow(muX_R, 2)+pow(muY_R, 2)+c1)*(varX_R+varY_R+c2));
    SSIM_G=((2*muX_G*muY_G+c1)*(2*varXY_G+c2))/((pow(muX_G, 2)+pow(muY_G, 2)+c1)*(varX_G+varY_G+c2));
    SSIM_B=((2*muX_B*muY_B+c1)*(2*varXY_B+c2))/((pow(muX_B, 2)+pow(muY_B, 2)+c1)*(varX_B+varY_B+c2));
    
    
    delete [] PixelArray;
    
    /////////////////////////OUTPUT STARTS HERE /////////////////////////
    
    //New image specs:
    cout<<"\nNew Image Specifications:"<<endl;
    cout << "Output Image Width="<<newWidth<<endl;
    cout << "Output Image Height="<<newHeight<<endl;
    cout << "\nPSNR Red="<<PSNR_R<<endl;
    cout << "PSNR Green="<<PSNR_G<<endl;
    cout << "PSNR Blue="<<PSNR_B<<endl;
    cout<<"PSNR Total="<<PSNR_R+PSNR_G+PSNR_B<<endl;
    
    cout << "\nMSE Red="<<MSE_R<<endl;
    cout << "MSE Green="<<MSE_G<<endl;
    cout << "MSE Blue="<<MSE_B<<endl;
    cout<<"MSE Total="<<MSE_R+MSE_G+MSE_B<<endl;
    cout << "\nPSNR Red="<<PSNR_R<<endl;
    cout << "PSNR Green="<<PSNR_G<<endl;
    cout << "PSNR Blue="<<PSNR_B<<endl;
    cout<<"PSNR Total="<<PSNR_R+PSNR_G+PSNR_B<<endl;
    cout << "\nSSIM Red="<<SSIM_R<<endl;
    cout << "SSIM Green="<<SSIM_G<<endl;
    cout << "SSIM Blue="<<SSIM_B<<endl;
    cout<<"SSIM Total="<<SSIM_R+SSIM_B+SSIM_G<<endl;
    
    
    
    // Write image data (filename specified by second argument) from image data matrix
    
    if (!(file=fopen(argv[2],"wb"))) {
        cout << "Cannot open file: " << argv[2] << endl;
        exit(1);
    }
    fwrite(tempImage, sizeof(unsigned char),  newHeight*newWidth*BytesPerPixel, file);
    fclose(file);
    return 0;
}
