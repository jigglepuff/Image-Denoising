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
    
    //Uncomment to debug: 
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
Function: LinearFilter
Input: Center pixel of interest (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), filter dimension (int) [eg. you would put 5 for a 5x5 filter]
Purpose: Applies linear/averaging filter on input image array to blur the image
Output: New intensity value for one pixel (int)
*/

int LinearFilter(Pixel center, Pixel* PixArray, int color, int height, int width, int windowDim)
{
    int average=0;
    int row =center.getRow();
    int col = center.getCol();
    int windowSize=windowDim*windowDim;
    int bounds=windowDim/2;
    vector<int> neighbors;
    if (row<bounds || col<bounds || row>height-1-bounds || col>width-1-bounds) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                average=(center.getRed())*windowSize;
                break;
            case 1:
                average=(center.getGreen())*windowSize;
                break;
            case 2:
                average=(center.getBlue())*windowSize;
                break;
                
            default:
                break;
        }
    }
    else//within bounds
    {
        neighbors=center.getNeighbors(width,windowDim);
        for (int i=0; i<windowSize; i++) {
            switch (color) {
                case 0:
                    average+=PixArray[neighbors[i]].getRed();
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
        for (int i=0; i<windowSize; i++) {
            neighbors.pop_back();
        }
    }
    return average/windowSize;
}


/*
Function: MedianFilter
Input: Center pixel of interest (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), filter dimension (int) [eg. you would put 5 for a 5x5 filter]
Purpose: Applies non-linear/median filter on input image array to blur the image
Output: New intensity value for one pixel (int)
*/

int MedianFilter(Pixel center, Pixel* PixArray, int color, int height, int width, int windowDim)
{
    int median=0;
    int row =center.getRow();
    int col = center.getCol();
    vector<int> neighbors;
    if (row==0 || col==0 || row==height-1 || col==width-1) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                median=(center.getRed());
                break;
            case 1:
                median=(center.getGreen());
                break;
            case 2:
                median=(center.getBlue());
                break;
            default:
                break;
        }
    }
    else//within bounds
    {
        neighbors=center.getNeighbors(width,windowDim); //Store pixel index of all neighbors in a vector
        for (int i=0; i<9; i++) { //use pixel index to get corresponding RGB intensity
            switch (color) {
                case 0:
                    neighbors.at(i)=PixArray[neighbors.at(i)].getRed();
                    break;
                case 1:
                    neighbors.at(i)=PixArray[neighbors.at(i)].getGreen();

                    break;
                    
                case 2:
                    neighbors.at(i)=PixArray[neighbors.at(i)].getBlue();
                    break;
                    
                default:
                    break;
            }
        }//end of index to intensity conversion
        
        sort(neighbors.begin(), neighbors.end());
        
        //Uncomment to debug:
        /*
        for (int i=0; i<9; i++) {
            cout<<i<<"="<<neighbors.at(i)<<" ";
        }
         */

        median=neighbors.at(4);
        
        //Uncomment to debug:
        //cout<<"\nmedian="<<median<<"\n"<<endl;

        for (int i=0; i<9; i++) {
            neighbors.pop_back();
        }
    }
    return median;
}


/*
Function: GaussianFilter()
Input: Center pixel of interest (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), Gaussian PDF variance (float), filter dimension (int) [eg. you would put 5 for a 5x5 filter]
Purpose: Applies Gaussian filter on input image array based on pixel intensities to blur the image
Output: New intensity value for one pixel (int)
*/
int GaussianFilter(Pixel center, Pixel* PixArray, int color, int height, int width, float variance, int windowDim)
{
    float weightedAvg=0;
    int row =center.getRow();
    int col = center.getCol();
    float gaussianRatio=0;
    float exponent=0;
    float gaussianSum=0;
    vector<int> neighbors;
    vector<float> weights;
    if (row==0 || col==0 || row==height-1 || col==width-1) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                weightedAvg=(center.getRed());
                break;
            case 1:
                weightedAvg=(center.getGreen());
                break;
            case 2:
                weightedAvg=(center.getBlue());
                break;
                
            default:
                break;
        }
    }
    else//within bounds
    {
        gaussianSum=0;
        exponent=0;
        gaussianRatio=0;
        for (int i=-1; i<2; i++) {
            for (int j=-1; j<2; j++) {
                exponent=-(pow(i,2)+pow(j,2))/(2*pow(variance, 2));
                gaussianRatio=(1/(2*3.141593*pow(variance, 2)))*exp(exponent);
                weights.push_back(gaussianRatio);//store calculated gaussian weights in vector
                gaussianSum+=gaussianRatio;  //Find sum of guassian weights for normalization
            }
        }
        
        neighbors=center.getNeighbors(width,windowDim); //Store pixel index of all neighbors in a vector
        weightedAvg=0;
        for (int i=0; i<9; i++) {
            switch (color) {
                case 0:
                    weightedAvg+=((weights.at(i))/gaussianSum)*float(PixArray[neighbors[i]].getRed());
                    break;
                case 1:
                    weightedAvg+=((weights.at(i))/gaussianSum)*float(PixArray[neighbors[i]].getGreen());
                    break;
                    
                case 2:
                    weightedAvg+=((weights.at(i))/gaussianSum)*float(PixArray[neighbors[i]].getBlue());
                    break;
                    
                default:
                    break;
            }
        }
    }
    
    for (int i=0; i<9; i++) {
        neighbors.pop_back();
        weights.pop_back();
    }
    
    weightedAvg=int(weightedAvg);
    return weightedAvg;
}

/*
Function: BilateralFilter()
Input: Center pixel of interest (Pixel), Image array (Pixel*), color R=0, G=1, B=2 (int), 
       height of image (int), width of image (int), variance of spatial gaussian function (float), variance of intensity gaussian function (float)
Purpose: Applies Gaussian filter on input image array based on surrounding pixels' spatial positions and intensities to blur the image
Output: New intensity value for one pixel (int)
*/
int BilateralFilter(Pixel center, Pixel* PixArray, int color, int height, int width, float spatialVar, float RGBVar)
{
    int row =center.getRow();
    int col = center.getCol();
    float weightedAvg=0;
    float GaussianSum=0;
    
    float GaussianRatio_spatial=0;
    float exponent_spatial=0;
    
    float GaussianRatio_RGB=0;
    float exponent_RGB=0;
    
    vector<int> neighbors;
    vector<float> weights_spatial;
    vector<float> weights_Total;
    
    if (row==0 || col==0 || row==height-1 || col==width-1) {//at least one neighbour is out of bounds
        switch (color) {
            case 0:
                weightedAvg=(center.getRed());
                break;
            case 1:
                weightedAvg=(center.getGreen());
                break;
            case 2:
                weightedAvg=(center.getBlue());
                break;
                
            default:
                break;
        }
    }
    else//within bounds
    {
        exponent_spatial=0;
        GaussianRatio_spatial=0;
        
        //Calculation of Spatial Weights
        for (int i=-1; i<2; i++) {
            for (int j=-1; j<2; j++) {
                exponent_spatial=-(pow(i,2)+pow(j,2))/(2*pow(spatialVar, 2));
                GaussianRatio_spatial=(1/(2*3.141593*pow(spatialVar, 2)))*exp(exponent_spatial);
                weights_spatial.push_back(GaussianRatio_spatial);//store calculated gaussian weights in vector
            }
        }//end of Calculation of Spatial Weights
        
        //Calculation of Intensity Weights
        GaussianSum=0;
        neighbors=center.getNeighbors(width,3); //Store pixel index of all neighbors in a vector
        for (int i=0; i<9; i++) {
            switch (color) {
                case 0:
                    exponent_RGB=-abs(PixArray[neighbors[i]].getRed()-center.getRed())/(2*pow(RGBVar,2));
                    GaussianRatio_RGB=(1/(2*3.141593*pow(RGBVar, 2)))*exp(exponent_RGB);
                    weights_Total.push_back(GaussianRatio_RGB*weights_spatial[i]);
                    GaussianSum+=GaussianRatio_RGB*weights_spatial[i];
                    break;
                case 1:
                    exponent_RGB=-abs(PixArray[neighbors[i]].getGreen()-center.getGreen())/(2*pow(RGBVar,2));
                    GaussianRatio_RGB=(1/(2*3.141593*pow(RGBVar, 2)))*exp(exponent_RGB);
                    weights_Total.push_back(GaussianRatio_RGB*weights_spatial[i]);
                    GaussianSum+=GaussianRatio_RGB*weights_spatial[i];
                    break;
                case 2:
                    exponent_RGB=-abs(PixArray[neighbors[i]].getBlue()-center.getBlue())/(2*pow(RGBVar,2));
                    GaussianRatio_RGB=(1/(2*3.141593*pow(RGBVar, 2)))*exp(exponent_RGB);
                    weights_Total.push_back(GaussianRatio_RGB*weights_spatial[i]);
                    GaussianSum+=GaussianRatio_RGB*weights_spatial[i];
                    break;
                default:
                    break;
            }
        }
        
        weightedAvg=0;
        for (int i=0; i<9; i++) {
            switch (color) {
                case 0:
                    weightedAvg+=(weights_Total.at(i)/GaussianSum)*float(PixArray[neighbors[i]].getRed());
                    break;
                case 1:
                    weightedAvg+=(weights_Total.at(i)/GaussianSum)*float(PixArray[neighbors[i]].getGreen());
                    break;
                    
                case 2:
                     weightedAvg+=(weights_Total.at(i)/GaussianSum)*float(PixArray[neighbors[i]].getBlue());
                    break;
                    
                default:
                    break;
            }
        }
    }
    

    for (int i=0; i<9; i++) {
        neighbors.pop_back();
        weights_spatial.pop_back();
        weights_Total.pop_back();
    }
    weightedAvg=int(weightedAvg);
    //Uncomment to debug:
    /*
    cout<<"red="<<center.getRed()<<endl;
    cout<<"Final weightedAvg="<<weightedAvg<<"\n"<<endl;
    */
    return weightedAvg;
}


/*
Function: main()
Input: command line input
*/
int main(int argc, const char * argv[]) {
    
   ///--------------------Command Line Inputs--------------//
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
    

    int newWidth;
    int newHeight;
    
    if (argc<6) {
        cout << "Syntax Error - Incorrect Parameter Usage:" << endl;
        cout << "program_name input_image.raw output_image.raw [BytesPerPixel] [InputImageWidth] [InputImageHeight]" << endl;
        exit(1);
    }

    ///-------------------- Output to Command Line------------------------//
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
    
   //-------------------- Set-up: Extract Image information ------------------------//
    //Initializing Output Image to all black pixels
    int pixelCount=0;
    Pixel* PixelArray;
    PixelArray = new Pixel [SizeWidth*SizeHeight];
    
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
    
    //-------------------- Uncomment filter functions you wish to implement------------------------//
    pixelCount=0;
    //PSNR variables
    float PSNR_R=0,PSNR_G=0,PSNR_B=0;
    float MSE_R=0,MSE_G=0,MSE_B=0;
    //SSIM variables
    float muX_R=0, muY_R=0,varX_R=0, varY_R=0, varXY_R=0;
    float muX_G=0, muY_G=0,varX_G=0, varY_G=0, varXY_G=0;
    float muX_B=0, muY_B=0,varX_B=0, varY_B=0, varXY_B=0;
   // filterType=atoi(argv[7]);
    
    for (int i=0; i<newHeight; i++)
    {
        for(int j=0; j<newWidth; j++)
        {
            for (int k=0; k<BytesPerPixel; k++) {
                
                //Comment & Uncomment the following lines to implement various filters: 

                /* Original Image */
                //tempImage[i][j][k]=originaldata[i][j][k];
                
                /* Unfiltered Image */
                //tempImage[i][j][k]=Imagedata[i][j][k];
               
                /* Mean 3x3 */
                tempImage[i][j][k]=char(computeNeighborhoodMean(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,3));
                
                /* Mean 7x7 */
                //tempImage[i][j][k]=char(computeNeighborhoodMean(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,7));
                
                /* Median 3x3 */
                //tempImage[i][j][k]=char(MedianFilter(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,3));
                
                /* Gaussian 3x3 */
                //tempImage[i][j][k]=char(GaussianFilter(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,1.6,3));
                
                /* Bilateral 3x3 */
                //tempImage[i][j][k]=char(BilateralFilter(PixelArray[pixelCount],PixelArray,k,SizeHeight,SizeWidth,1.6,15));
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
    
    //PSNR Calculation
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
    
    //Output to command line:
    
    //New image specs:
    cout<<"\nNew Image Specifications:"<<endl;
    cout << "Output Image Width="<<newWidth<<endl;
    cout << "Output Image Height="<<newHeight<<endl;

    cout<<"MSE Total="<<MSE_R+MSE_G+MSE_B<<endl;
    cout<<"PSNR Total="<<PSNR_R+PSNR_G+PSNR_B<<endl;
    cout<<"SSIM Total="<<SSIM_R+SSIM_G+SSIM_B<<endl;

    // Write image data (filename specified by second argument) from image data matrix
    
    if (!(file=fopen(argv[2],"wb"))) {
        cout << "Cannot open file: " << argv[2] << endl;
        exit(1);
    }
    fwrite(tempImage, sizeof(unsigned char),  newHeight*newWidth*BytesPerPixel, file);
    fclose(file);
    return 0;
}
