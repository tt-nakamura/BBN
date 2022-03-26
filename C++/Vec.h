// W. H. Press, et al, "Numerical Recipes"

#ifndef __Vec_h__
#define __Vec_h__

template <class T>
class Vec {
private:
    int nn;	// size of array. upper index is nn-1
    T *v;
public:
    Vec();
    explicit Vec(int n);		// Zero-based array
    Vec(const T &a, int n);	//initialize to constant value
    Vec(const T *a, int n);	// Initialize to array
    Vec(const Vec &rhs);	// Copy constructor
    Vec & operator=(const Vec &rhs);	//assignment
    Vec & operator=(const T &a);	//assign a to every element
    inline T & operator[](const int i);	//i'th element
    inline const T & operator[](const int i) const;
    inline int size() const;
    ~Vec();
};

template <class T>
Vec<T>::Vec() : nn(0), v(0) {}

template <class T>
Vec<T>::Vec(int n) : nn(n), v(new T[n]) {}

template <class T>
Vec<T>::Vec(const T& a, int n) : nn(n), v(new T[n])
{
    for(int i=0; i<n; i++)
        v[i] = a;
}

template <class T>
Vec<T>::Vec(const T *a, int n) : nn(n), v(new T[n])
{
    for(int i=0; i<n; i++)
        v[i] = *a++;
}

template <class T>
Vec<T>::Vec(const Vec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
    for(int i=0; i<nn; i++)
        v[i] = rhs[i];
}

template <class T>
Vec<T> & Vec<T>::operator=(const Vec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if (this != &rhs)
    {
        if (nn != rhs.nn) {
            if (v != 0) delete [] (v);
            nn=rhs.nn;
            v= new T[nn];
        }
        for (int i=0; i<nn; i++)
            v[i]=rhs[i];
    }
    return *this;
}

template <class T>
Vec<T> & Vec<T>::operator=(const T &a)	//assign a to every element
{
    for (int i=0; i<nn; i++)
        v[i]=a;
    return *this;
}

template <class T>
inline T & Vec<T>::operator[](const int i)	//subscripting
{
    return v[i];
}

template <class T>
inline const T & Vec<T>::operator[](const int i) const	//subscripting
{
    return v[i];
}

template <class T>
inline int Vec<T>::size() const
{
    return nn;
}

template <class T>
Vec<T>::~Vec()
{
    if (v != 0)
        delete[] (v);
}

typedef Vec<double> Vec_DP, Vec_O_DP, Vec_IO_DP;
typedef Vec<int> Vec_INT, Vec_O_INT, Vec_IO_INT;
typedef const Vec<double> Vec_I_DP;
typedef const Vec<int> Vec_I_INT;

#endif // __Vec_h__