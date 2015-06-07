class CalcError 
{ 
private: 
	int nx;
	double *dm; 
	double dmax; 
public: 
	explicit CalcError(double *tdm, int tnx): 
						dm(tdm), nx(tnx), dmax(0) {} 
	CalcError(const CalcError& m, split): 
						dm(m.dm), nx(m.nx), dmax(m.dmax) {} 
	void operator()(const blocked_range<int>& r) 
	{ 
		int begin = r.begin(), end = r.end(); 
		for ( int i=begin; i<end; i++ ) 
			{ 
				if ( dmax < dm[i] ) dmax = dm[i]; 
				
				//printf("%i: dmax=%f\n",i,dmax);
			}
	} 
	void join(const CalcError& m) 
	{ 
		if ( dmax < m.dmax ) dmax = m.dmax;
	} 
	double Result() { return dmax; }
};
