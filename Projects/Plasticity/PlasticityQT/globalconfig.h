#ifndef GLOBALCONFIG_H
#define GLOBALCONFIG_H

#include <QDockWidget>

namespace Ui {
class GlobalConfig;
}

class GlobalConfig : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit GlobalConfig(QWidget *parent = 0);
    ~GlobalConfig();
    
private:
    Ui::GlobalConfig *ui;
};

#endif // GLOBALCONFIG_H
