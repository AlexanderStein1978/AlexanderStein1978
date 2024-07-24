

class SoundVector;

class SoundMatrix
{
public:
    SoundMatrix(const int numColumns, const int numRows);
    ~SoundMatrix();

    SoundVector& operator[](const int column);
    SoundVector operator*(const SoundVector& right);

private:
    SoundVector* mData;
    int mColumnCount, mRowCount;
};
