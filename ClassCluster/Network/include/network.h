#ifndef NETWORK_H
#define NETWORK_H


#include <QTimer>
#include <QMap>

class QTcpSocket;

class Window;


class Network : public QObject
{
    Q_OBJECT
public:
    enum Command
    {
        START, STOP, SEND_FLAGS, RESET, MOVE, TRIGGER_SNAP_SHOT, WRITE_SNAP_SHOT, RESTORE_SNAP_SHOT, GET_SETTINGS_AND_POTENTIALS, SET_SETTINGS, SET_POTENTIAL, RELOAD_POTENTIALS, STOP_CALC,
        ROTATE, DATA_RECEIVED, DATA_FOLLOWING, LOGMESSAGE, ERROR_INCOMPLETE, ERROR_UNKNOWN_COMMAND, NO_KNOWN_COMMAND
    };

    Network(Window *window);
    virtual ~Network();
    void SendCommand(const QByteArray& command);
    void SendCommand(const Command command);
    void SendCommand(const Command command, const QByteArray& data);
    void SendFlags(const char flags);

protected:
    virtual void dataReceived();
    virtual void commandReceived(const Command command);
    double readDouble(bool complete);
    quint32 readUint32(bool complete);
    QString readString();

    static const int SIZE_OF_COMMAND_STRINGS;
    static const qint64 SIZE_OF_LOG_MESSAGE_TIME;

    qint64 minimumDataToRead;
    QTcpSocket* mSocket;
    Window* mWindow;
    QMap<QByteArray, Command> mCommandMap;
    QTimer mTimer;

    std::vector<QTcpSocket*> mOldSockets;

private slots:
    void run();

private:
    void RecoverFromError(char* receivedBuffer);
    bool ReadStreamedObject(QByteArray& buffer);

    int mResendCount;
    QByteArray mLastSentCommand;
};

#endif // NETWORK_H
