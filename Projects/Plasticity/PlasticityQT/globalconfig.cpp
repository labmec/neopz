#include "globalconfig.h"
#include "ui_globalconfig.h"

GlobalConfig::GlobalConfig(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::GlobalConfig)
{
    ui->setupUi(this);
}

GlobalConfig::~GlobalConfig()
{
    delete ui;
}
