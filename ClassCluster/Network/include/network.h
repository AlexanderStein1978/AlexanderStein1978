#ifndef NETWORK_H
#define NETWORK_H


#include <QThread>
#include <QMap>

class QTcpSocket;

class Window;


class Network : public QThread
{
public:
    enum Command
    {
        START, STOP, RESET, MOVE, TRIGGER_SNAP_SHOT, WRITE_SNAP_SHOT, RESTORE_SNAP_SHOT, SET_KINETIC_ENERGY, SET_POTENTIAL_RANGE_SCALE, SET_SPEED, SET_STEP_SIZE,
        RELOAD_POTENTIALS, STOP_CALC, ROTATE, SET_LAYER_DISTANCE, DATA_RECEIVED, DATA_FOLLOWING, ERROR_INCOMPLETE, ERROR_UNKNOWN_COMMAND
    };

    Network(Window *window);
    virtual ~Network();
    void SendCommand(const QByteArray& command);
    void SendCommand(const Command command);

protected:
    void run();
    virtual void dataReceived();
    virtual void commandReceived(const Command command) = 0;
    double readDouble(bool complete);
    quint32 readUint32(bool complete);

    qint64 minimumDataToRead;
    QTcpSocket* mSocket;
    Window* mWindow;
    QMap<QByteArray, Command> mCommandMap;

    static const int SIZE_OF_COMMAND_STRINGS;

private:
    bool continueRunning;
};

#endif // NETWORK_H
