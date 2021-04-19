

class VariableSet
{
public:
	virtual ~VariableSet() { }
	virtual void print(IloCplex &cplex) = 0;
};

class MtzVariableSet : public VariableSet
{
public:
	~MtzVariableSet();
	void print(IloCplex &cplex);

	IloBoolVarArray x;
	IloBoolVarArray v;
	IloIntVarArray u;
};