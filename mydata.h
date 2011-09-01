typedef struct pathNode
{
	Point pos;
	int val;
	
	int m_iFscore;
	int m_iHscore;
	int m_iGscore;
	
	pathNode* pParent;
	
	pathNode(int x = 0,int y = 0)
	{
		pos.row = x;
		pos.col= y;
		m_iFscore = 0;
		m_iGscore = 0;
		m_iBscore = 0;
		pParent = NULL;
		//val = 0;
	}
}
pathNode,*Nodeptr;