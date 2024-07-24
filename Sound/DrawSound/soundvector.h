
class SoundVector
{
public:
    SoundVector();
    SoundVector(const int numElements);
    SoundVector(const SoundVector& temp);
    ~SoundVector();

    void resize(const int newSize);
    double& operator[](const int index);
    const double& operator[](const int index) const;
    int getSize() const;
    SoundVector sigmoid() const;
    int getIndexOfMax() const;

private:
    double *mData;
    int mSize;
};
