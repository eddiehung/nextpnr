/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  Miodrag Milanovic <miodrag@symbioticeda.com>
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include <QAction>
#include <QFileDialog>
#include <QGridLayout>
#include <QIcon>
#include <QSplitter>
#include "designwidget.h"
#include "fpgaviewwidget.h"
#include "jsonparse.h"
#include "log.h"
#include "mainwindow.h"
#include "ScintillaEdit.h"
#include "SciLexer.h"
#include "pythontab.h"

static void initBasenameResource() { Q_INIT_RESOURCE(base); }

NEXTPNR_NAMESPACE_BEGIN

BaseMainWindow::BaseMainWindow(QWidget *parent) : QMainWindow(parent), ctx(nullptr)
{
    initBasenameResource();
    qRegisterMetaType<std::string>();

    log_files.clear();
    log_streams.clear();

    setObjectName(QStringLiteral("BaseMainWindow"));
    resize(1024, 768);

    createMenusAndBars();

    QWidget *centralWidget = new QWidget(this);

    QGridLayout *gridLayout = new QGridLayout(centralWidget);
    gridLayout->setSpacing(6);
    gridLayout->setContentsMargins(11, 11, 11, 11);

    QSplitter *splitter_h = new QSplitter(Qt::Horizontal, centralWidget);
    QSplitter *splitter_v = new QSplitter(Qt::Vertical, splitter_h);
    splitter_h->addWidget(splitter_v);

    gridLayout->addWidget(splitter_h, 0, 0, 1, 1);

    setCentralWidget(centralWidget);

    DesignWidget *designview = new DesignWidget();
    designview->setMinimumWidth(300);
    designview->setMaximumWidth(300);
    splitter_h->addWidget(designview);

    connect(this, SIGNAL(contextChanged(Context *)), designview, SLOT(newContext(Context *)));

    connect(designview, SIGNAL(info(std::string)), this, SLOT(writeInfo(std::string)));

    tabWidget = new QTabWidget();
#ifndef NO_PYTHON
    PythonTab *pythontab = new PythonTab();
    tabWidget->addTab(pythontab, "Python");
    connect(this, SIGNAL(contextChanged(Context *)), pythontab, SLOT(newContext(Context *)));
    connect(this, SIGNAL(executePython(QString)), pythontab, SIGNAL(execute(QString)));
    connect(this, SIGNAL(runPythonScript(QString)), pythontab, SIGNAL(runScript(QString)));
#endif
    info = new InfoTab();
    tabWidget->addTab(info, "Info");

    centralTabWidget = new QTabWidget();
    centralTabWidget->setTabsClosable(true);

    FPGAViewWidget *fpgaView = new FPGAViewWidget();

    connect(centralTabWidget, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));

    centralTabWidget->addTab(fpgaView, "Graphics");
    centralTabWidget->tabBar()->tabButton(0, QTabBar::RightSide)->resize(0, 0);

    connect(this, SIGNAL(contextChanged(Context *)), fpgaView, SLOT(newContext(Context *)));

    splitter_v->addWidget(centralTabWidget);
    splitter_v->addWidget(tabWidget);
}

BaseMainWindow::~BaseMainWindow() {}

void BaseMainWindow::closeTab(int index)
{
    centralTabWidget->removeTab(index);
}

void BaseMainWindow::writeInfo(std::string text) { info->info(text); }

void BaseMainWindow::createMenusAndBars()
{
    actionNew = new QAction("New", this);
    QIcon iconNew;
    iconNew.addFile(QStringLiteral(":/icons/resources/new.png"));
    actionNew->setIcon(iconNew);
    actionNew->setShortcuts(QKeySequence::New);
    actionNew->setStatusTip("New project file");
    connect(actionNew, SIGNAL(triggered()), this, SLOT(new_proj()));

    actionOpen = new QAction("Open", this);
    QIcon iconOpen;
    iconOpen.addFile(QStringLiteral(":/icons/resources/open.png"));
    actionOpen->setIcon(iconOpen);
    actionOpen->setShortcuts(QKeySequence::Open);
    actionOpen->setStatusTip("Open an existing project file");
    connect(actionOpen, SIGNAL(triggered()), this, SLOT(open_proj()));

    QAction *actionSave = new QAction("Save", this);
    QIcon iconSave;
    iconSave.addFile(QStringLiteral(":/icons/resources/save.png"));
    actionSave->setIcon(iconSave);
    actionSave->setShortcuts(QKeySequence::Save);
    actionSave->setStatusTip("Save existing project to disk");
    connect(actionSave, SIGNAL(triggered()), this, SLOT(save_proj()));
    actionSave->setEnabled(false);

    QAction *actionExit = new QAction("Exit", this);
    QIcon iconExit;
    iconExit.addFile(QStringLiteral(":/icons/resources/exit.png"));
    actionExit->setIcon(iconExit);
    actionExit->setShortcuts(QKeySequence::Quit);
    actionExit->setStatusTip("Exit the application");
    connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));

    QAction *actionAbout = new QAction("About", this);

    menuBar = new QMenuBar();
    menuBar->setGeometry(QRect(0, 0, 1024, 27));
    QMenu *menu_File = new QMenu("&File", menuBar);
    QMenu *menu_Help = new QMenu("&Help", menuBar);
    menuBar->addAction(menu_File->menuAction());
    menuBar->addAction(menu_Help->menuAction());
    setMenuBar(menuBar);

    mainToolBar = new QToolBar();
    addToolBar(Qt::TopToolBarArea, mainToolBar);

    documentsToolBar = new QToolBar();
    addToolBar(Qt::TopToolBarArea, documentsToolBar);

    statusBar = new QStatusBar();
    setStatusBar(statusBar);

    menu_File->addAction(actionNew);
    menu_File->addAction(actionOpen);
    menu_File->addAction(actionSave);
    menu_File->addSeparator();
    menu_File->addAction(actionExit);
    menu_Help->addAction(actionAbout);

    mainToolBar->addAction(actionNew);
    mainToolBar->addAction(actionOpen);
    mainToolBar->addAction(actionSave);

    QAction *actionDocNew = new QAction("New Document", this);
    QIcon iconDocNew;
    iconDocNew.addFile(QStringLiteral(":/icons/resources/page.png"));
    actionDocNew->setIcon(iconDocNew);
    actionDocNew->setStatusTip("New document");
    connect(actionDocNew, SIGNAL(triggered()), this, SLOT(new_doc()));

    QAction *actionDocOpen = new QAction("Open Document", this);
    QIcon iconDocOpen;
    iconDocOpen.addFile(QStringLiteral(":/icons/resources/page_edit.png"));
    actionDocOpen->setIcon(iconDocOpen);
    actionDocOpen->setStatusTip("Open document");

    QAction *actionDocSave = new QAction("Open Document", this);
    QIcon iconDocSave;
    iconDocSave.addFile(QStringLiteral(":/icons/resources/page_save.png"));
    actionDocSave->setIcon(iconDocSave);
    actionDocSave->setStatusTip("Save document");

    QAction *actionDocExecute = new QAction("Execute", this);
    QIcon iconDocExecute;
    iconDocExecute.addFile(QStringLiteral(":/icons/resources/page_go.png"));
    actionDocExecute->setIcon(iconDocExecute);
    actionDocExecute->setStatusTip("Execute document");    
    connect(actionDocExecute, SIGNAL(triggered()), this, SLOT(execute_doc()));

    QAction *actionDocRun = new QAction("Execute", this);
    QIcon iconDocRun;
    iconDocRun.addFile(QStringLiteral(":/icons/resources/page_lightning.png"));
    actionDocRun->setIcon(iconDocRun);
    actionDocRun->setStatusTip("Run from file");    
    connect(actionDocRun, SIGNAL(triggered()), this, SLOT(run_doc()));
#ifdef NO_PYTHON
    actionDocRun->setEnabled(false);
#endif

    documentsToolBar->addAction(actionDocNew);
    documentsToolBar->addAction(actionDocOpen);
    documentsToolBar->addAction(actionDocSave);
    documentsToolBar->addAction(actionDocExecute);
    documentsToolBar->addAction(actionDocRun);
}

void BaseMainWindow::new_doc()
{
    ScintillaEdit *editor = new ScintillaEdit();
    editor->setLexer(SCLEX_PYTHON);
    editor->styleClearAll();
	editor->setMarginWidthN(0, 35);
	editor->setScrollWidth(200);
	editor->setScrollWidthTracking(1);

    editor->styleSetFore(SCE_P_DEFAULT, 0x000000);
    editor->styleSetFore(SCE_P_COMMENTLINE, 0x008000);
    editor->styleSetFore(SCE_P_NUMBER, 0xFF0000);
    editor->styleSetFore(SCE_P_STRING, 0x808080);
    editor->styleSetFore(SCE_P_CHARACTER, 0x808080);
    editor->styleSetFore(SCE_P_WORD, 0x0000FF);
    editor->styleSetFore(SCE_P_TRIPLE, 0xFF8000);
    editor->styleSetFore(SCE_P_TRIPLEDOUBLE, 0xFF8000);
    editor->styleSetFore(SCE_P_CLASSNAME, 0x000000);
    editor->styleSetFore(SCE_P_DEFNAME, 0xFF00FF);
    editor->styleSetFore(SCE_P_OPERATOR, 0x000080);
    editor->styleSetFore(SCE_P_IDENTIFIER, 0x000000);
    editor->styleSetFore(SCE_P_COMMENTBLOCK, 0x008000);
	editor->styleSetFore(SCE_P_DECORATOR, 0xFF8000);
    
    centralTabWidget->addTab(editor, "New");
    centralTabWidget->setCurrentIndex(centralTabWidget->count()-1);
}

void BaseMainWindow::execute_doc()
{
#ifndef NO_PYTHON
   ScintillaEdit* editor = static_cast<ScintillaEdit *>(centralTabWidget->currentWidget());
   QString data = QString(editor->get_doc()->get_char_range(0,editor->get_doc()->length()));
   Q_EMIT executePython(data);
#endif   
}

void BaseMainWindow::run_doc()
{
#ifndef NO_PYTHON
    QString fileName = QFileDialog::getOpenFileName(this, QString("Run Python"), QString(), QString("*.py"));
    if (!fileName.isEmpty()) {
        tabWidget->setCurrentIndex(0);
        Q_EMIT runPythonScript(fileName);
    }
#endif   
}
NEXTPNR_NAMESPACE_END
