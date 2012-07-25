/**
 *@brief bvh_parser     : parse bvh data
         dataTransform  : transform bvh data(joint rotation angles,offset) to 3-dimensional points,
         3D points can be used to draw the motion

         recommend you rewrite the code with one-dimensional array
 *@author songtianyi630@163.com
 */


#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>

/////////////////////////////////////////////////////////////////
/**
 *stack ,C++, zero-based
 *clear() push() top() empty() pop() size()
 *free()
 */
typedef int T;
struct myStack{
    int curr,size_limit;
    T *array;
    myStack(int s){/*input stack size*/
        curr = -1; size_limit = s;
        array = new T[s];
    }
    void clear(){curr = -1;}/*empty stack*/
    void free(){ if(array != 0) delete [] array; array = 0; clear();}/*free stack*/
    T top(){
        if( curr > -1 && curr < size_limit ){return array[curr];}
        else return -1;}
    bool push(T v){
        if(curr + 1 < size_limit){ array[++curr] = v; return true;}
        return false;}
    bool empty(){
        if(curr < 0)return true;
        else return false;}
    void pop(){curr--;}
    int  size(){return curr + 1;}
    ~myStack()
    {
        assert(array == 0);/*release memory explicitly*/
    }
};
//////////////////////////////////////////////////////////////////////////


///define
#define mymin(a,b) (a>b?b:a)
#define mymax(a,b) (a>b?a:b)
#define MAX_JOINT 100
#define MAX_FRAME 1000
#define MAX_COL (3 + MAX_JOINT*3)

///global varibles
struct joint
{
    int id;//one-based,[1,MAX_JOINT]
    int channel;//channel {3,6}
    int channels[6];//{0,1,2} - {X,Y,Z}
    char joint_name[10];
    double offset[3];//offset value

}joint_hiry[ MAX_JOINT + 1 ];


//bvh parser
int joint_num;//joint count
int joint_end;//END SITE count
int channels_num;//total channels
int frame_num;//frame count
double frame_time;//time per frame

int parent_of[MAX_JOINT+1];//parent of joint, parent_of[j] = -1 denote j is illegal, = 0 denote j is root joint
bool is_end[MAX_JOINT + 1];//mark end site joint
double mat[ MAX_FRAME ][ MAX_COL ];//original MOTION data


void bvh_parser(const char *bvh_dir)
{
    FILE *bvh = fopen(bvh_dir,"r");
    if(bvh == NULL)
    {
        printf("[open file error] %s\n",bvh_dir);
        return;
    }

    //declarition
    char buffer[256];
    myStack my(MAX_JOINT+1);

    //initialization
    joint_num = 0;
    joint_end = 0;
    channels_num = 0;
    memset(parent_of,-1,sizeof(parent_of));
    memset(is_end,0,sizeof(is_end));
    parent_of[1] = 0;

    //hierarchy construct
    while( fscanf(bvh,"%s",buffer) != EOF )
    {
        if( strcmp(buffer,"{") == 0 )
        {
            my.push(joint_num);
        }
        else if( strcmp(buffer,"}") == 0 )
        {
            int c = my.top(); my.pop();

            if(c == 1)
            {
                //back to root
                break;
            }

            parent_of[ c ] = my.top();

        }
        else if( strcmp(buffer,"OFFSET") == 0 )
        {
            fscanf(bvh,"%lf%lf%lf",&joint_hiry[joint_num].offset[0],\
                   &joint_hiry[joint_num].offset[1],\
                   &joint_hiry[joint_num].offset[2]);
        }
        else if( strcmp(buffer,"CHANNELS") == 0)
        {
            fscanf(bvh,"%d",&joint_hiry[joint_num].channel);
            channels_num += joint_hiry[joint_num].channel;
            for(int i = 0;i < joint_hiry[joint_num].channel; i++)
            {
                fscanf(bvh,"%s",buffer);
                //x = 0, y = 1, z = 2
                //position rotation | rotation
                if( buffer[0] == 'X')
                {
                    joint_hiry[joint_num].channels[i] = 0;
                }
                else if( buffer[0] == 'Y' )
                {
                    joint_hiry[joint_num].channels[i] = 1;
                }
                else if( buffer[0] == 'Z')
                {
                    joint_hiry[joint_num].channels[i] = 2;
                }
            }
        }
        else if( strcmp(buffer,"JOINT") == 0 || strcmp(buffer,"ROOT") == 0 || strcmp(buffer,"End") == 0)
        {
            fscanf(bvh,"%s",joint_hiry[++joint_num].joint_name);

            if(buffer[0] == 'E')
            {
                is_end[joint_num] = true;
                joint_end++;
            }
        }
    }
    my.free();
    //end of construction


    fscanf(bvh,"%s",buffer);fscanf(bvh,"%s",buffer);
    fscanf(bvh,"%d",&frame_num);
    fscanf(bvh,"%s",buffer);fscanf(bvh,"%s",buffer);
    fscanf(bvh,"%lf",&frame_time);

    int col = 6 + (joint_num - joint_end - 1)*3;
    assert(col == channels_num);
    assert(frame_num < MAX_FRAME);
    assert(joint_num < MAX_JOINT);
    for(int i = 0;i < frame_num;i++)
    {
        for(int j = 0;j < channels_num;j++)
        {
            fscanf(bvh,"%lf",&mat[i][j]);
        }
    }

    fclose( bvh ); bvh = 0;
}


///define
#define PI 3.14159265358979323846

///
const int DOUBLE_SIZE = sizeof(double);
const double unitM[4][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };//unit matrix
double M[ MAX_JOINT + 1][4][4];//transformation matrix

double Tmat[ MAX_FRAME ][MAX_JOINT+1][ 3 ];//transformed data(3D points)



/**
 *matrix multiplication
 * a[4][4] X b[4][4] = rs[4][4]
 * a[4][4] X b[4][1] = rs[4][1]
 *rs will be overwrite
 */
void matrixMulti444(double rs[][4],double a[][4],double b[][4])
{
    double temp[4][4]; memset(temp,0,sizeof(temp));
    for(int i = 0;i < 4;i++)
    {
        for(int j = 0;j < 4;j++)
        {
            for(int k = 0;k < 4;k++)
            {
                temp[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    memcpy(rs,temp,sizeof(temp));
}
void matrixMulti441(double rs[],double a[][4],double b[])
{
    double temp[4]; memset(temp,0,sizeof(temp));
    for(int i = 0;i < 4;i++)
    {
        for(int k = 0;k < 4;k++)
        {
            temp[i] += a[i][k]*b[k];
        }
    }
    memcpy(rs,temp,sizeof(temp));
}



void getTransMatrix(double rs[][4],int f,int j,int end_num)
{
    ///T*R

    //rs = T
    double T[4][4] = {
        {1,0,0,joint_hiry[j].offset[0]},
        {0,1,0,joint_hiry[j].offset[1]},
        {0,0,1,joint_hiry[j].offset[2]},
        {0,0,0,1}};
    memcpy(rs,T,sizeof(T));

    //rs = rs*R
    int i = 0;
    if(joint_hiry[j].channel == 6)
    {
        //root joint
        for(;i < 3;i++)
        {
            //position
            if(joint_hiry[j].channels[i] == 0)
            {
                //Xposition
                double t[4][4] = {
                    {1,0,0,mat[f][i]},
                    {0,1,0,0},
                    {0,0,1,0},
                    {0,0,0,1} };
                matrixMulti444(rs,rs,t);
            }
            else if(joint_hiry[j].channels[i] == 1)
            {
                //Yposition
                double t[4][4] = {
                    {1,0,0,0},
                    {0,1,0,mat[f][i]},
                    {0,0,1,0},
                    {0,0,0,1} };
                matrixMulti444(rs,rs,t);
            }
            else{
                //Zposition
                double t[4][4] = {
                    {1,0,0,0},
                    {0,1,0,0},
                    {0,0,1,mat[f][i]},
                    {0,0,0,1} };
                matrixMulti444(rs,rs,t);
            }
        }
    }

    int off_col = (j-end_num)*3;
    for(int k = 0;k < 3;i++,k++)
    {
        //rotation
        double x = mat[f][ off_col + k]*PI/180;
        if(joint_hiry[j].channels[i] == 0)
        {
            //Xrotation
            double t[4][4] = { {1,0,0,0},{0,cos(x),-sin(x),0},{0,sin(x),cos(x),0},{0,0,0,1} };
            matrixMulti444(rs,rs,t);
        }
        else if(joint_hiry[j].channels[i] == 1)
        {
            //Yrotation
            double t[4][4] = { {cos(x),0,sin(x),0},{0,1,0,0},{-sin(x),0,cos(x),0},{0,0,0,1} };
            matrixMulti444(rs,rs,t);
        }
        else{
            //Zrotation
            double t[4][4] = { {cos(x),-sin(x),0,0},{sin(x),cos(x),0,0},{0,0,1,0},{0,0,0,1} };
            matrixMulti444(rs,rs,t);
        }
    }
}


void dataTransform()
{
    memcpy(M[0],unitM,sizeof(unitM));
    for(int i = 0;i < frame_num;i++)
    {
        int end_num = 0;
        for(int j = 1;j <= joint_num;j++)
        {
            if(is_end[j] == false)
            {
                //M[parent_of[j]]*R(channels[0])*R()...
                double trans[4][4];
                getTransMatrix(trans,i,j,end_num);
                matrixMulti444(M[j], M[parent_of[j]], trans);
            }
            else end_num++;

            //S[4][1] = {0,0,0,1}
            //E[4][1] = {offset[0],offset[1],offset[2],1};
            //E = M[parent_of[j]]*E  S = M[parent_of[j]]*S
            if(j == 1)
            {
                double S[4] = {0,0,0,1};
                matrixMulti441(S,M[1],S);
                memcpy(Tmat[i][0],S,3*DOUBLE_SIZE);
            }
            else{

                double E[4] = {joint_hiry[j].offset[0],joint_hiry[j].offset[1],joint_hiry[j].offset[2],1};
                matrixMulti441(E,M[ parent_of[j] ],E);
                memcpy(Tmat[i][j-1],E,3*DOUBLE_SIZE);
            }
        }
    }
}
int main(){return 0;}
