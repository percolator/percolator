#ifndef VIRTUALSET_H_
#define VIRTUALSET_H_


class VirtualSet : public DataSet
{
public:
	VirtualSet(DataSet &mother,int fold,int ix);
	virtual ~VirtualSet();
};

#endif /*VIRTUALSET_H_*/
