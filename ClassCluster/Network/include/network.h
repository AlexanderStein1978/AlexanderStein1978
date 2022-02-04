#ifndef NETWORK_H
#define NETWORK_H


#include <QThread>

class QTcpSocket;

class Window;


class Network : private QThread
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
    virtual void dataReceived() = 0;

    QTcpSocket* mSocket;
    Window* mWindow;
    QMap<QByteArray, Command> mCommandMap;

    static const int SIZE_OF_COMMAND_STRINGS;

private:
    bool continueRunning;
};

#endif // NETWORK_H
