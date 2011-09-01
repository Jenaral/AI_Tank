#include "Tank.h"
#include <string.h>
#include <stdio.h>
//请勿修改以上头文件
/*

您可以在这里添加您所需头文件

*/
//#include <windows.h>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include <utility>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#define Z_NORMAL 0
#define Z_STAR 1
#define Z_PRIME 2
//matrix.h
template <class T>
class Matrix {
public:
	Matrix();
	Matrix(int rows, int columns);
	Matrix(const Matrix<T> &other);
	Matrix<T> & operator= (const Matrix<T> &other);
	~Matrix();
	// all operations except product modify the matrix in-place.
	void resize(int rows, int columns);
	void identity(void);
	void clear(void);
	T& operator () (int x, int y);
	T trace(void);
	Matrix<T>& transpose(void);
	Matrix<T> product(Matrix<T> &other);
	int minsize(void) {
		return ((m_rows < m_columns) ? m_rows : m_columns);
	}
	int columns(void) {
		return m_columns;
	}
	int rows(void) {
		return m_rows;
	}
private:
	T **m_matrix;
	int m_rows;
	int m_columns;
};
//end matrix.h

//munkres.h


class Munkres {
public:
	void solve(Matrix<double> &m);
private:
	inline bool find_uncovered_in_matrix(double,int&,int&);
	inline bool pair_in_list(const std::pair<int,int> &, const std::list<std::pair<int,int> > &);
	int step1(void);
	int step2(void);
	int step3(void);
	int step4(void);
	int step5(void);
	int step6(void);
	Matrix<int> mask_matrix;
	Matrix<double> matrix;
	bool *row_mask;
	bool *col_mask;
	int saverow, savecol;
};
//end munkres.h
//matrix.cpp


/*export*/ template <class T>
Matrix<T>::Matrix() {
	m_rows = 0;
	m_columns = 0;
	m_matrix = NULL;
}

/*export*/ template <class T>
Matrix<T>::Matrix(const Matrix<T> &other) {
	if ( other.m_matrix != NULL ) {
		// copy arrays
		m_matrix = NULL;
		resize(other.m_rows, other.m_columns);
		for ( int i = 0 ; i < m_rows ; i++ )
			for ( int j = 0 ; j < m_columns ; j++ )
				m_matrix[i][j] = other.m_matrix[i][j];
	} else {
		m_matrix = NULL;
		m_rows = 0;
		m_columns = 0;
	}
}

/*export*/ template <class T>
Matrix<T>::Matrix(int rows, int columns) {
	m_matrix = NULL;
	resize(rows, columns);
}

/*export*/ template <class T>
Matrix<T> &
Matrix<T>::operator= (const Matrix<T> &other) {
	if ( other.m_matrix != NULL ) {
		// copy arrays
		resize(other.m_rows, other.m_columns);
		for ( int i = 0 ; i < m_rows ; i++ )
			for ( int j = 0 ; j < m_columns ; j++ )
				m_matrix[i][j] = other.m_matrix[i][j];
	} else {
		// free arrays
		for ( int i = 0 ; i < m_columns ; i++ )
			delete [] m_matrix[i];

		delete [] m_matrix;

		m_matrix = NULL;
		m_rows = 0;
		m_columns = 0;
	}

	return *this;
}

/*export*/ template <class T>
Matrix<T>::~Matrix() {
	if ( m_matrix != NULL ) {
		// free arrays
		for ( int i = 0 ; i < m_rows ; i++ )
			delete [] m_matrix[i];

		delete [] m_matrix;
	}
	m_matrix = NULL;
}

/*export*/ template <class T>
void
Matrix<T>::resize(int rows, int columns) {
	if ( m_matrix == NULL ) {
		// alloc arrays
		m_matrix = new T*[rows]; // rows
		for ( int i = 0 ; i < rows ; i++ )
			m_matrix[i] = new T[columns]; // columns

		m_rows = rows;
		m_columns = columns;
		clear();
	} else {
		// save array pointer
		T **new_matrix;
		// alloc new arrays
		new_matrix = new T*[rows]; // rows
		for ( int i = 0 ; i < rows ; i++ ) {
			new_matrix[i] = new T[columns]; // columns
			for ( int j = 0 ; j < columns ; j++ )
				new_matrix[i][j] = 0;
		}

		// copy data from saved pointer to new arrays
		int minrows = std::min<int>(rows, m_rows);
		int mincols = std::min<int>(columns, m_columns);
		for ( int x = 0 ; x < minrows ; x++ )
			for ( int y = 0 ; y < mincols ; y++ )
				new_matrix[x][y] = m_matrix[x][y];

		// delete old arrays
		if ( m_matrix != NULL ) {
			for ( int i = 0 ; i < m_rows ; i++ )
				delete [] m_matrix[i];

			delete [] m_matrix;
		}

		m_matrix = new_matrix;
	}

	m_rows = rows;
	m_columns = columns;
}

/*export*/ template <class T>
void
Matrix<T>::identity() {
	assert( m_matrix != NULL );

	clear();

	int x = std::min<int>(m_rows, m_columns);
	for ( int i = 0 ; i < x ; i++ )
		m_matrix[i][i] = 1;
}

/*export*/ template <class T>
void
Matrix<T>::clear() {
	assert( m_matrix != NULL );

	for ( int i = 0 ; i < m_rows ; i++ )
		for ( int j = 0 ; j < m_columns ; j++ )
			m_matrix[i][j] = 0;
}

/*export*/ template <class T>
T 
Matrix<T>::trace() {
	assert( m_matrix != NULL );

	T value = 0;

	int x = std::min<int>(m_rows, m_columns);
	for ( int i = 0 ; i < x ; i++ )
		value += m_matrix[i][i];

	return value;
}

/*export*/ template <class T>
Matrix<T>& 
Matrix<T>::transpose() {
	assert( m_rows > 0 );
	assert( m_columns > 0 );

	int new_rows = m_columns;
	int new_columns = m_rows;

	if ( m_rows != m_columns ) {
		// expand matrix
		int m = std::max<int>(m_rows, m_columns);
		resize(m,m);
	}

	for ( int i = 0 ; i < m_rows ; i++ ) {
		for ( int j = i+1 ; j < m_columns ; j++ ) {
			T tmp = m_matrix[i][j];
			m_matrix[i][j] = m_matrix[j][i];
			m_matrix[j][i] = tmp;
		}
	}

	if ( new_columns != new_rows ) {
		// trim off excess.
		resize(new_rows, new_columns);
	}

	return *this;
}

/*export*/ template <class T>
Matrix<T> 
Matrix<T>::product(Matrix<T> &other) {
	assert( m_matrix != NULL );
	assert( other.m_matrix != NULL );
	assert ( m_columns == other.m_rows );

	Matrix<T> out(m_rows, other.m_columns);

	for ( int i = 0 ; i < out.m_rows ; i++ ) {
		for ( int j = 0 ; j < out.m_columns ; j++ ) {
			for ( int x = 0 ; x < m_columns ; x++ ) {
				out(i,j) += m_matrix[i][x] * other.m_matrix[x][j];
			}
		}
	}

	return out;
}

/*export*/ template <class T>
T&
Matrix<T>::operator ()(int x, int y) {
	assert ( x >= 0 );
	assert ( y >= 0 );
	assert ( x < m_rows );
	assert ( y < m_columns );
	assert ( m_matrix != NULL );
	return m_matrix[x][y];
}

//end matrix.cpp
//munkres.cpp


bool 
Munkres::find_uncovered_in_matrix(double item, int &row, int &col) {
	for ( row = 0 ; row < matrix.rows() ; row++ )
		if ( !row_mask[row] )
			for ( col = 0 ; col < matrix.columns() ; col++ )
				if ( !col_mask[col] )
					if ( matrix(row,col) == item )
						return true;

	return false;
}

bool 
Munkres::pair_in_list(const std::pair<int,int> &needle, const std::list<std::pair<int,int> > &haystack) {
	for ( std::list<std::pair<int,int> >::const_iterator i = haystack.begin() ; i != haystack.end() ; i++ ) {
		if ( needle == *i )
			return true;
	}

	return false;
}

int 
Munkres::step1(void) {
	for ( int row = 0 ; row < matrix.rows() ; row++ )
		for ( int col = 0 ; col < matrix.columns() ; col++ )
			if ( matrix(row,col) == 0 ) {
				bool isstarred = false;
				for ( int nrow = 0 ; nrow < matrix.rows() ; nrow++ )
					if ( mask_matrix(nrow,col) == Z_STAR )
						isstarred = true;

				if ( !isstarred ) {
					for ( int ncol = 0 ; ncol < matrix.columns() ; ncol++ )
						if ( mask_matrix(row,ncol) == Z_STAR )
							isstarred = true;
				}

				if ( !isstarred ) {
					mask_matrix(row,col) = Z_STAR;
				}
			}

			return 2;
}

int 
Munkres::step2(void) {
	int covercount = 0;
	for ( int row = 0 ; row < matrix.rows() ; row++ )
		for ( int col = 0 ; col < matrix.columns() ; col++ )
			if ( mask_matrix(row,col) == Z_STAR ) {
				col_mask[col] = true;
				covercount++;
			}

			int k = matrix.minsize();

			if ( covercount >= k ) {
#ifdef DEBUG
				std::cout << "Final cover count: " << covercount << std::endl;
#endif
				return 0;
			}

#ifdef DEBUG
			std::cout << "Munkres matrix has " << covercount << " of " << k << " Columns covered:" << std::endl;
			for ( int row = 0 ; row < matrix.rows() ; row++ ) {
				for ( int col = 0 ; col < matrix.columns() ; col++ ) {
					std::cout.width(8);
					std::cout << matrix(row,col) << ",";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
#endif


			return 3;
}

int 
Munkres::step3(void) {
	/*
	Main Zero Search

	1. Find an uncovered Z in the distance matrix and prime it. If no such zero exists, go to Step 5
	2. If No Z* exists in the row of the Z', go to Step 4.
	3. If a Z* exists, cover this row and uncover the column of the Z*. Return to Step 3.1 to find a new Z
	*/
	if ( find_uncovered_in_matrix(0, saverow, savecol) ) {
		mask_matrix(saverow,savecol) = Z_PRIME; // prime it.
	} else {
		return 5;
	}

	for ( int ncol = 0 ; ncol < matrix.columns() ; ncol++ )
		if ( mask_matrix(saverow,ncol) == Z_STAR ) {
			row_mask[saverow] = true; //cover this row and
			col_mask[ncol] = false; // uncover the column containing the starred zero
			return 3; // repeat
		}

		return 4; // no starred zero in the row containing this primed zero
}

int 
Munkres::step4(void) {
	std::list<std::pair<int,int> > seq;
	// use saverow, savecol from step 3.
	std::pair<int,int> z0(saverow, savecol);
	std::pair<int,int> z1(-1,-1);
	std::pair<int,int> z2n(-1,-1);
	seq.insert(seq.end(), z0);
	int row, col = savecol;
	/*
	Increment Set of Starred Zeros

	1. Construct the ``alternating sequence'' of primed and starred zeros:

	Z0 : Unpaired Z' from Step 4.2 
	Z1 : The Z* in the column of Z0
	Z[2N] : The Z' in the row of Z[2N-1], if such a zero exists 
	Z[2N+1] : The Z* in the column of Z[2N]

	The sequence eventually terminates with an unpaired Z' = Z[2N] for some N.
	*/
	bool madepair;
	do {
		madepair = false;
		for ( row = 0 ; row < matrix.rows() ; row++ )
			if ( mask_matrix(row,col) == Z_STAR ) {
				z1.first = row;
				z1.second = col;
				if ( pair_in_list(z1, seq) )
					continue;

				madepair = true;
				seq.insert(seq.end(), z1);
				break;
			}

			if ( !madepair )
				break;

			madepair = false;

			for ( col = 0 ; col < matrix.columns() ; col++ )
				if ( mask_matrix(row,col) == Z_PRIME ) {
					z2n.first = row;
					z2n.second = col;
					if ( pair_in_list(z2n, seq) )
						continue;
					madepair = true;
					seq.insert(seq.end(), z2n);
					break;
				}
	} while ( madepair );

	for ( std::list<std::pair<int,int> >::iterator i = seq.begin() ;
		i != seq.end() ;
		i++ ) {
			// 2. Unstar each starred zero of the sequence.
			if ( mask_matrix(i->first,i->second) == Z_STAR )
				mask_matrix(i->first,i->second) = Z_NORMAL;

			// 3. Star each primed zero of the sequence,
			// thus increasing the number of starred zeros by one.
			if ( mask_matrix(i->first,i->second) == Z_PRIME )
				mask_matrix(i->first,i->second) = Z_STAR;
	}

	// 4. Erase all primes, uncover all columns and rows, 
	for ( int row = 0 ; row < mask_matrix.rows() ; row++ )
		for ( int col = 0 ; col < mask_matrix.columns() ; col++ )
			if ( mask_matrix(row,col) == Z_PRIME )
				mask_matrix(row,col) = Z_NORMAL;

	for ( int i = 0 ; i < matrix.rows() ; i++ ) {
		row_mask[i] = false;
	}

	for ( int i = 0 ; i < matrix.columns() ; i++ ) {
		col_mask[i] = false;
	}

	// and return to Step 2. 
	return 2;
}

int 
Munkres::step5(void) {
	/*
	New Zero Manufactures

	1. Let h be the smallest uncovered entry in the (modified) distance matrix.
	2. Add h to all covered rows.
	3. Subtract h from all uncovered columns
	4. Return to Step 3, without altering stars, primes, or covers. 
	*/
	double h = 0;
	for ( int row = 0 ; row < matrix.rows() ; row++ ) {
		if ( !row_mask[row] ) {
			for ( int col = 0 ; col < matrix.columns() ; col++ ) {
				if ( !col_mask[col] ) {
					if ( (h > matrix(row,col) && matrix(row,col) != 0) || h == 0 ) {
						h = matrix(row,col);
					}
				}
			}
		}
	}

	for ( int row = 0 ; row < matrix.rows() ; row++ )
		for ( int col = 0 ; col < matrix.columns() ; col++ ) {
			if ( row_mask[row] )
				matrix(row,col) += h;

			if ( !col_mask[col] )
				matrix(row,col) -= h;
		}

		return 3;
}

void 
Munkres::solve(Matrix<double> &m) {
	// Linear assignment problem solution
	// [modifies matrix in-place.]
	// matrix(row,col): row major format assumed.

	// Assignments are remaining 0 values
	// (extra 0 values are replaced with -1)
#ifdef DEBUG
	std::cout << "Munkres input matrix:" << std::endl;
	for ( int row = 0 ; row < m.rows() ; row++ ) {
		for ( int col = 0 ; col < m.columns() ; col++ ) {
			std::cout.width(8);
			std::cout << m(row,col) << ",";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
#endif

	bool notdone = true;
	int step = 1;

	this->matrix = m;
	// Z_STAR == 1 == starred, Z_PRIME == 2 == primed
	mask_matrix.resize(matrix.rows(), matrix.columns());

	row_mask = new bool[matrix.rows()];
	col_mask = new bool[matrix.columns()];
	for ( int i = 0 ; i < matrix.rows() ; i++ ) {
		row_mask[i] = false;
	}

	for ( int i = 0 ; i < matrix.columns() ; i++ ) {
		col_mask[i] = false;
	}

	while ( notdone ) {
		switch ( step ) {
			case 0:
				notdone = false;
				break;
			case 1:
				step = step1();
				break;
			case 2:
				step = step2();
				break;
			case 3:
				step = step3();
				break;
			case 4:
				step = step4();
				break;
			case 5:
				step = step5();
				break;
		}
	}

	// Store results
	for ( int row = 0 ; row < matrix.rows() ; row++ )
		for ( int col = 0 ; col < matrix.columns() ; col++ )
			if ( mask_matrix(row,col) == Z_STAR )
				matrix(row,col) = 0;
			else
				matrix(row,col) = -1;

#ifdef DEBUG
	std::cout << "Munkres output matrix:" << std::endl;
	for ( int row = 0 ; row < matrix.rows() ; row++ ) {
		for ( int col = 0 ; col < matrix.columns() ; col++ ) {
			std::cout.width(1);
			std::cout << matrix(row,col) << ",";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
#endif

	m = matrix;

	delete [] row_mask;
	delete [] col_mask;
}
//end munkres.cpp
/*

您可以在这里添加您的自定义函数

*/

short pToi[MAP_HEIGHT+1][MAP_WIDTH+1];
Point iTop[MAP_HEIGHT*MAP_WIDTH+1];
short dist[MAP_HEIGHT*MAP_WIDTH+1][MAP_HEIGHT*MAP_WIDTH+1];
Point next[10];//the next step for each tank
bool lastMove[10];
Point base;
std::vector<Point> nextStep[10];
short dir[4][2]={0,-1,0,1,-1,0,1,0};//up down left right
OrderType direction[4]={GOLEFT,GORIGHT,GOUP,GODOWN};
Point nearestS[10];//the NeutralSource assinged for each tank
//char* ds[4]={"left","right","up","down"};
std::vector<Point> bS,rS,nS;
typedef enum {SHOOT,DRIVE}STATE;
STATE tankState[10];
DataForAI mydata;


bool isInSet(Point& x,std::vector<Point>& openSet)
{
	for(std::vector<Point>::iterator i=openSet.begin();i!=openSet.end();i++)
	{
		if(i->row==x.row&&i->col==x.col)
			return true;
	}
	return false;
}

short heuristic_estimate_of_distance(Point& s,Point& t)
{
	//return dist[pToi[s.row][s.col]][pToi[t.row][t.col]];
	return abs(s.row-t.row)+abs(s.col-t.col);

}
bool isNotCollision(Point& t,DataForAI& data)
{
	for(int i=0;i<data.myID;i++)
	{
		if(t.row==next[i].row&&t.col==next[i].col)
			return false;
	}
	//if(data.map[t.row][t.col].whoIsHere!=-1)
	//return false;
	return true;
}
Point aStar(Point& start,Point& goal,DataForAI& data)
{

	std::vector<Point> closeSet,openSet;
	//if(data.round==1&&data.myID==4)
	//printf("astar for tank%d:start(%d,%d),goal(%d,%d)\n",data.myID,start.row,start.col,goal.row,goal.col);
	Point cameFrom[MAP_HEIGHT+2][MAP_WIDTH+2];
	short g_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	short h_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	short f_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	Point origin;
	origin.row=-1;
	origin.col=-1;
	if(goal.row<0||goal.col<0)
		return origin;
	for(short i=0;i<MAP_HEIGHT+2;i++)
		for(short j=0;j<MAP_WIDTH+2;j++)
		{
			cameFrom[i][j]=origin;//no path
			g_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
			h_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
			f_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
		}
		openSet.push_back(start);//set the open set containing the initial node
		g_score[start.row][start.col]=0;
		h_score[start.row][start.col]=heuristic_estimate_of_distance(start,goal);
		f_score[start.row][start.col]=h_score[start.row][start.col];

		while(!openSet.empty())
		{
			short min;
			Point x;
			std::vector<Point>::iterator iter=openSet.begin();
			x.row=iter->row;
			x.col=iter->col;
			//get the node which have the lowest fscore in open set
			for(min=f_score[iter->row][iter->col];iter!=openSet.end();iter++)
			{
				if(min>=f_score[iter->row][iter->col])
				{
					min=f_score[iter->row][iter->col];
					x.row=iter->row;
					x.col=iter->col;
				}

			}
			//if(data.round==1&&data.myID==4)
			//printf("round%d tank%d:x(%d,%d) have lowest f_score:%d\n",data.round,data.myID,x.row,x.col,f_score[x.row][x.col]);
			if(x.row==goal.row&&x.col==goal.col)
			{
				//printf("next step:%d,%d\n",cameFrom[goal.row][goal.col].row,cameFrom[goal.row][goal.col].col);
				return cameFrom[goal.row][goal.col];
			}
			//return reconstructPath(cameFrom[goal.row][goal.col]);
			//remove x from open set
			//if(data.round==1&&data.myID==4)
			//printf("round%d tank%d:erase x(%d,%d) from open set and add to close set\n",data.round,data.myID,x.row,x.col);
			for(std::vector<Point>::iterator i=openSet.begin();i!=openSet.end();i++)
			{
				if(x.row==i->row&&x.col==i->col)
				{
					openSet.erase(i);
					break;
				}
			}
			//add x to close set
			closeSet.push_back(x);
			short r=x.row,c=x.col;//index of x
			//each neighbor y of x
			for(short j=0;j<4;j++)
			{
				//printf("direction %d:%d,%d\n",j,dir[j][0],dir[j][1]);
				short nr=r+dir[j][0],nc=c+dir[j][1];
				Point y;//index of y
				y.row=nr;
				y.col=nc;			 

				if(nr<=MAP_HEIGHT && nr>=1 && nc<=MAP_WIDTH && nc>=1)//bound detection
					if(data.map[nr][nc].type!=STONE&&isNotCollision(y,data))//not stone
					{
						short mt=data.map[nr][nc].whoIsHere;
						if(mt==-1||mt<=data.myID)
						{
							if(isInSet(y,closeSet))
								continue;
							short tentative_g_score=g_score[r][c]+1;//tentative_g_score := g_score[x] + dist_between(x,y)
							bool tentative_is_better;
							if(!isInSet(y,openSet))
							{
								openSet.push_back(y);
								tentative_is_better=true;

							}
							else if(tentative_g_score<g_score[nr][nc])
								tentative_is_better=true;
							else
								tentative_is_better=false;
							if(tentative_is_better==true)
							{
								/*
								came_from[y] := x
								g_score[y] := tentative_g_score
								h_score[y] := heuristic_estimate_of_distance(y, goal)
								f_score[y] := g_score[y] + h_score[y]
								*/


								cameFrom[nr][nc]=x;
								if(!isInSet(x,nextStep[data.myID]))
									nextStep[data.myID].push_back(x);

								//printf("(%d,%d)<-(%d,%d)\n",nr,nc,x.row,x.col);
								g_score[nr][nc]=tentative_g_score;
								h_score[nr][nc]=heuristic_estimate_of_distance(y,goal);
								f_score[nr][nc]=g_score[nr][nc]+h_score[nr][nc];
							}
						}//if obstacles
					}
			}//end for 4 directions
		}//end while
		//if(data.round==1)
		//printf("end astar for tank%d\n",data.myID);
		return origin;
}//end astar algorithm
bool isNeighbor(short s,short t)
{
	short d=abs(s-t);
	if(d==MAP_HEIGHT||d==1)
		return true;
	return false;
}

void floyd(DataForAI& data)
{
	//initialize the graph
	for(int i=1;i<=MAP_HEIGHT*MAP_WIDTH;i++)
		for(int j=1;j<=MAP_HEIGHT*MAP_WIDTH;j++)
		{
			if(isNeighbor(i,j))
			{//is neibors
				if(data.map[iTop[i].row][iTop[i].col].type==STONE||data.map[iTop[j].row][iTop[j].col].type==STONE)
					dist[i][j]=MAP_HEIGHT*MAP_WIDTH;
				else if(data.map[iTop[i].row][iTop[i].col].type==BRICK||data.map[iTop[j].row][iTop[j].col].type==BRICK)
					dist[i][j]=2;
				else if(data.map[iTop[i].row][iTop[i].col].type==BREAKBRICK||data.map[iTop[j].row][iTop[j].col].type==BREAKBRICK)
					dist[i][j]=1;
				else
					dist[i][j]=1;
			}
			else
				dist[i][j]=MAP_HEIGHT*MAP_WIDTH;
		}



		//loop of floyd algorithm
		for(int k=1;k<=MAP_HEIGHT*MAP_WIDTH;k++)
			for(int i=1;i<=MAP_HEIGHT*MAP_WIDTH;i++)
				for(int j=1;j<=MAP_HEIGHT*MAP_WIDTH;j++)
					if(dist[i][k]+dist[k][j]<dist[i][j])
						dist[i][j]=dist[i][k]+dist[k][j];
		//printf("floyd complete!\n");
		/*
		for(int i=1;i<=MAP_HEIGHT*MAP_WIDTH;i++)
		{
		for(int j=1;j<=MAP_HEIGHT*MAP_WIDTH;j++)
		printf("d(%d,%d)=%d ",i,j,dist[i][j]);
		printf("\n");
		}
		*/

}
short estAstarDistance(Point& start,Point& goal,DataForAI& data)
{
	std::vector<Point> closeSet,openSet,step;
	short g_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	short h_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	short f_score[MAP_HEIGHT+2][MAP_WIDTH+2];
	for(short i=0;i<MAP_HEIGHT+2;i++)
		for(short j=0;j<MAP_WIDTH+2;j++)
		{
			//cameFrom[i][j]=origin;//no path
			g_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
			h_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
			f_score[i][j]=MAP_HEIGHT*MAP_WIDTH;
		}
		openSet.push_back(start);//set the open set containing the initial node
		g_score[start.row][start.col]=0;
		h_score[start.row][start.col]=heuristic_estimate_of_distance(start,goal);
		f_score[start.row][start.col]=h_score[start.row][start.col];

		while(!openSet.empty())
		{
			short min;
			Point x;
			std::vector<Point>::iterator iter=openSet.begin();
			x.row=iter->row;
			x.col=iter->col;
			//get the node which have the lowest fscore in open set
			for(min=f_score[iter->row][iter->col];iter!=openSet.end();iter++)
			{
				if(min>=f_score[iter->row][iter->col])
				{
					min=f_score[iter->row][iter->col];
					x.row=iter->row;
					x.col=iter->col;
				}

			}
			//if(data.round==1&&data.myID==4)
			//printf("round%d tank%d:x(%d,%d) have lowest f_score:%d\n",data.round,data.myID,x.row,x.col,f_score[x.row][x.col]);
			if(x.row==goal.row&&x.col==goal.col)
				return step.size();
			//return reconstructPath(cameFrom[goal.row][goal.col]);
			//remove x from open set
			//if(data.round==1&&data.myID==4)
			//printf("round%d tank%d:erase x(%d,%d) from open set and add to close set\n",data.round,data.myID,x.row,x.col);
			for(std::vector<Point>::iterator i=openSet.begin();i!=openSet.end();i++)
			{
				if(x.row==i->row&&x.col==i->col)
				{
					openSet.erase(i);
					break;
				}
			}
			//add x to close set
			closeSet.push_back(x);
			short r=x.row,c=x.col;//index of x
			//each neighbor y of x
			for(short j=0;j<4;j++)
			{
				//printf("direction %d:%d,%d\n",j,dir[j][0],dir[j][1]);
				short nr=r+dir[j][0],nc=c+dir[j][1];
				Point y;//index of y
				y.row=nr;
				y.col=nc;			 

				if(nr<=MAP_HEIGHT && nr>=1 && nc<=MAP_WIDTH && nc>=1)//bound detection
					if(data.map[nr][nc].type!=STONE)//not stone
					{
						if(isInSet(y,closeSet))
							continue;
						short tentative_g_score=g_score[r][c]+1;//tentative_g_score := g_score[x] + dist_between(x,y)
						bool tentative_is_better;
						if(!isInSet(y,openSet))
						{
							openSet.push_back(y);
							tentative_is_better=true;

						}
						else if(tentative_g_score<g_score[nr][nc])
							tentative_is_better=true;
						else
							tentative_is_better=false;
						if(tentative_is_better==true)
						{
							/*
							came_from[y] := x
							g_score[y] := tentative_g_score
							h_score[y] := heuristic_estimate_of_distance(y, goal)
							f_score[y] := g_score[y] + h_score[y]
							*/

							if(!isInSet(x,step))
								step.push_back(x);
							g_score[nr][nc]=tentative_g_score;
							h_score[nr][nc]=heuristic_estimate_of_distance(y,goal);
							f_score[nr][nc]=g_score[nr][nc]+h_score[nr][nc];
						}
					}//if obstacles
			}//end for 4 directions
		}//end while
		//if(data.round==1)
		//printf("end astar for tank%d\n",data.myID);
		return MAP_HEIGHT*MAP_WIDTH;

}
short estDistance(const Point& s,const Point& t)
{
	//return dist[pToi[s.row][s.col]][pToi[t.row][t.col]];
	return abs(s.row-t.row)+abs(s.col-t.col);


}
short estMyPower(Point& t,DataForAI& data)
{
		short sum=0;
	for(int i=-SOURCE_SIGHT;i<=SOURCE_SIGHT;i++)
	{
		for(int j=-SOURCE_SIGHT;j<=SOURCE_SIGHT;j++)
		{
			if(abs(i)+abs(j)<=SOURCE_SIGHT)
			{
				short nr=t.row+i,nc=t.col+j;
				if(data.myFlag=RED)
				{
					short enemy=data.map[nr][nc].whoIsHere;
					if(enemy>=0&&enemy<5)
						sum+=data.tank[enemy].attack+data.tank[enemy].life+data.tank[enemy].range;
				}
				else if(data.myFlag=BLUE)
				{
					short enemy=data.map[nr][nc].whoIsHere;
					if(enemy>5)
						sum+=data.tank[enemy].attack+data.tank[enemy].life+data.tank[enemy].range;
				}
			}

		}
	}
	return sum;
}

short estDangerous(const Point& t,DataForAI& data)
{
	short sum=0;
	for(int i=-SOURCE_SIGHT;i<=SOURCE_SIGHT;i++)
	{
		for(int j=-SOURCE_SIGHT;j<=SOURCE_SIGHT;j++)
		{
			if(abs(i)+abs(j)<=SOURCE_SIGHT)
			{
				short nr=t.row+i,nc=t.col+j;
				if(data.myFlag=RED)
				{
					short enemy=data.map[nr][nc].whoIsHere;
					if(enemy>5)
						sum+=data.tank[enemy].attack+data.tank[enemy].life+data.tank[enemy].range;
				}
				else if(data.myFlag=BLUE)
				{
					short enemy=data.map[nr][nc].whoIsHere;
					if(enemy>=0&&enemy<=4)
						sum+=data.tank[enemy].attack+data.tank[enemy].life+data.tank[enemy].range;
				}
			}

		}
	}
	return sum;
}
short estDistance(short a1,short a2,short b1,short b2)
{
	short sum=abs(a1-b1)+abs(a2-b2);
	return sum;

}
bool isInMap(short a,short b)
{
	if(a>0||a<MAP_HEIGHT||b>0||b<MAP_WIDTH)
		return true;
	return false;
}
void initPointIndex()
{
	for(short r=1;r<=MAP_HEIGHT;r++)
		for(short c=1;c<=MAP_WIDTH;c++)
		{
			pToi[r][c]=(r-1)*MAP_HEIGHT+c;
			//printf("%d,",(r-1)*MAP_HEIGHT+c);
		}

		for(short i=1;i<=MAP_HEIGHT*MAP_WIDTH;i++)
		{

			iTop[i].row=i/MAP_HEIGHT+1;
			iTop[i].col=i%MAP_WIDTH;
		}
}
void initState()
{
	for(int i=0;i<10;i++)
		tankState[i]==DRIVE;
}
bool isObstacle(Point& a,std::vector<Point>& nStep)
{
	std::vector<Point>::iterator i=nStep.begin();
	for(;i!=nStep.end();i++)
		if((a.row==i->row)&&(a.col==i->col))
			return true;
	return false;
}
void updateResource(DataForAI& data)
{
	bS.resize(0);
	rS.resize(0);
	nS.resize(0);
	for(short i=0;i<data.totalSource;i++)
	{
		short r=data.source[i].row;
		short c=data.source[i].col;
		if(data.map[r][c].isHereSource==1)//red
			rS.push_back(data.source[i]);
		if(data.map[r][c].isHereSource==2)//blue
			bS.push_back(data.source[i]);
		if(data.map[r][c].isHereSource==3)//belongs to nobody
			nS.push_back(data.source[i]);
	}
	std::vector<Point>::iterator iter;
	if(data.myFlag==RED)
		for(iter=bS.begin();iter!=bS.end();iter++)
			nS.push_back(*iter);
	else
		for(iter=rS.begin();iter!=rS.end();iter++)
			nS.push_back(*iter);



}
bool nearToTank(const Point& a,const Point& b)
{
	short as=0,bs=0;
	if(mydata.myFlag=RED)
		for(int i=0;i<5;i++)
		{
			as+=estDistance(a.row,a.col,mydata.tank[i].row,mydata.tank[i].col);
			bs+=estDistance(b.row,b.col,mydata.tank[i].row,mydata.tank[i].col);
		}
	else if(mydata.myFlag=BLUE)
		for(int i=5;i<10;i++)
		{
			as+=estDistance(a.row,a.col,mydata.tank[i].row,mydata.tank[i].col);
			bs+=estDistance(b.row,b.col,mydata.tank[i].row,mydata.tank[i].col);
		}

		if(as<bs)
			return true;
		return false;
}
void assignResource(DataForAI& data)
{

	//updateResource(data);
	std::vector<Point>::iterator iter=nS.begin();
	sort(nS.begin(),nS.end(),nearToTank);
	short m=nS.size(),n=5;
	short k=m<n?m:n;
	Point me;
	me.row=data.tank[data.myID].row;
	me.col=data.tank[data.myID].col;
	Matrix<double> kmMatrix(k,k);
	//printf("initialize the matrix\n");
	//initialize the matrix
	//initKmMatrix(kmMatrix);
	int i=0;
	for(iter=nS.begin();i<k;i++,iter++)//soure id
	{
		for(int j=0;j<k;j++)//tank id
		{

			kmMatrix(i,j)=1/double(estAstarDistance(me,*iter,data));
		}
		//i++;
	}
	//printf("apply the munkres matching \n");
	//apply the munkres matching 
	Munkres km;
	km.solve(kmMatrix);

	//printf(" Display solved matrix\n");
#ifdef DISP_KMRESULT
	// Display solved matrix.
	if(data.round==1)
	{
		for ( int row = 0 ; row < k ; row++ ) {
			for ( int col = 0 ; col < k ; col++ ) {
				std::cout.width(2);
				std::cout << kmMatrix(row,col) << ",";
			}
			std::cout << std::endl;
		}
	}
#endif
	//std::vector<Point>::iterator 
	iter=nS.begin();
	for (int r=0;r<k;r++) 
	{
		for (int c=0;c<k;c++) 
		{
			if(kmMatrix(r,c)==0)
			{
				nearestS[r]=*(iter+c);
				//if(data.round==1)
				printf("nS(%d,%d) assigned to tank%d\n",(*(iter+c)).row,(*(iter+c)).col,r);

			}
		}

	}
}
bool isAssigned(Point& s,DataForAI& data)
{
	short i;
	if(data.myFlag==RED)
	{
		for(i=0;i<5;i++)
			if(s.row==nearestS[i].row&&s.col==nearestS[i].col)
				return true;
	}
	else if(data.myFlag==BLUE)
	{
		for(i=5;i<10;i++)
			if(s.row==nearestS[i].row&&s.col==nearestS[i].col)
				return true;
	}

	return false;
}

void getNearestRS(DataForAI& data)
{
	updateResource(data);
	short min=MAP_HEIGHT+MAP_WIDTH;
	Point me;
	me.row=data.tank[data.myID].row;
	me.col=data.tank[data.myID].col;
	std::vector<Point>::iterator p;
	//assert(nS.empty());
	for(std::vector<Point>::iterator i=nS.begin();i!=nS.end();i++)
	{
		printf("(%d,%d),",i->row,i->col);
	}
	printf("\n");
	for(std::vector<Point>::iterator i=nS.begin();i!=nS.end();i++)
	{
		short tmp=estDistance(me,*i);//+estDangerous(*i,data)-estMyPower(me,data);
		if(min>tmp&&!isAssigned(*i,data))
		{
			nearestS[data.myID]=*i;		
		}
	}
}
void getNext(DataForAI& data)
{
	nextStep[data.myID].resize(0);
	short i=data.myID;
	Point goalTank;
	goalTank.row=data.tank[i].row;
	goalTank.col=data.tank[i].col;
	next[i]=aStar(nearestS[i],goalTank,data);
	//printf("round%d next step for tank%d:(%d,%d)\n",data.round,data.myID,next[i].row,next[i].col);
}

//平台0回合时调用此函数获取AI名称及坦克类型信息，请勿修改此函数声明。
extern "C" InitiateInfo chooseType()
{
	InitiateInfo Info;
	Info.tank[0]=Sniper;
	Info.tank[1]=Sniper;
	Info.tank[2]=Sniper;
	Info.tank[3]=Sniper;
	Info.tank[4]=Pioneer;

	strcpy(Info.aiName,"MyAIName"); //AI名请勿使用中文。 
	return Info;
}


//平台从第1回合开始调用此函数获得每回合指令，请勿修改此函数声明。
extern "C" Order makeOrder(DataForAI data)
{
	//printf("===================round%d tank%d==================\n",data.round,data.myID);
	Order order;
	order.type=STOP;
	mydata=data;
	updateResource(data);
	if(data.round==1&&data.myID==0)
	{
		assignResource(data);
		//getNext(data);
	}

	short range=data.tank[data.myID].range;
	//printf("range:%d\n",range);
	short vision=range;
	Point me;
	me.row=data.tank[data.myID].row;
	me.col=data.tank[data.myID].col;
	for(int i=-vision;i<=vision;i++)
	{
		for(int j=-vision;j<=vision;j++)
		{
			if(abs(i)+abs(j)<=vision)//in vision
			{
				short nr=me.row+i,nc=me.col+j;
				Point target;
				target.row=nr;
				target.col=nc;
				if(isInMap(nr,nc))//in map
					if(data.map[nr][nc].whoIsHere!=-1)//tank detected
					{
						short n=data.map[nr][nc].whoIsHere;
						if(data.tank[n].flag!=data.myFlag)
						{
							//printf("enemy tank%d\n",n);
							order.row=nr;
							order.col=nc;
							order.type=FIRE;
							return order;
						}

					}
					else if((data.map[nr][nc].type==BREAKBRICK||data.map[nr][nc].type==BRICK)&&isObstacle(target,nextStep[data.myID]))//bricks detected
					{
						//printf("tank%dbricks\n",data.myID);
						order.row=nr;
						order.col=nc;
						order.type=FIRE;
						return order;

					}
			}
		}
	}

	//

	//printf("the nearest point for tank%d (%d,%d)\n",data.myID,nearestS[data.myID].row,nearestS[data.myID].col);
	//assert(nearestS[data.myID].row>0);
	if(data.tank[data.myID].life!=0)
	{
		clock_t tstart,tend;
		float timeuse;
		tstart=clock();

		if(data.map[nearestS[data.myID].row][nearestS[data.myID].col].isHereSource==RedSource||data.map[data.tank[data.myID].row][data.tank[data.myID].col].isHereSource==RedSource||nextStep[data.myID].empty())
		{
			getNearestRS(data);
		}
		getNext(data);
		tend=clock();
		timeuse=(tend-tstart)/CLOCKS_PER_SEC;
		//printf("round%d tank%dchange resource(%d,%d) and next step:(%d,%d),%f seconds used\n",data.round,data.myID,nearestS[data.myID].row,nearestS[data.myID].col,next[data.myID].row,next[data.myID].col,timeuse);
	};

	/*
	if(data.map[next[data.myID].row][next[data.myID].col].whoIsHere!=-1||next[data.myID].row==-1)
	{
	while(1)
	{
	srand((int)time(0));
	short i=rand()%4;
	short nr=data.tank[data.myID].row+dir[i][0],nc=data.tank[data.myID].col+dir[i][1];
	if(data.map[nr][nc].type==PERVIOUS&&data.map[nr][nc].whoIsHere!=-1)
	{
	next[data.myID].row=nr;
	next[data.myID].col=nc;
	break;
	}
	}
	}
	*/
	short x=next[data.myID].row-data.tank[data.myID].row,y=next[data.myID].col-data.tank[data.myID].col;
	//printf("(x,y) for tank%d\n",x,y,data.myID);
	/*if(x<0)x=-1;
	else if(x>0)x=1;
	if(y<0)y=-1;
	else if(y>0)y=1;*/
	for(int i=0;i<4;i++)
	{
		if(dir[i][0]==x&&dir[i][1]==y)
		{

			order.type=direction[i];
			//printf("tank%d is going %s\n",data.myID,ds[i]);
			return order;
		}
		else
		{	
			order.type=STOP;
		}
		//random step
		/*
		else
		{
		srand((int)time(0));
		order.type=direction[rand()%4];
		return order;

		}*/
	}
	return order;
}


