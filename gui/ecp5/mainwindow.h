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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "../basewindow.h"

NEXTPNR_NAMESPACE_BEGIN

class MainWindow : public BaseMainWindow
{
    Q_OBJECT

  public:
    explicit MainWindow(std::unique_ptr<Context> context, ArchArgs args, QWidget *parent = 0);
    virtual ~MainWindow();

  public:
    void createMenu();
    void load_base_config(std::string filename);

  protected:
    void onDisableActions() override;
    void onJsonLoaded() override;
    void onRouteFinished() override;
    void onProjectLoaded() override;

  protected Q_SLOTS:
    void new_proj() override;
    void newContext(Context *ctx);
    void open_lpf();
    void open_base();
    void save_config();

  private:
    QAction *actionLoadLPF;
    QAction *actionLoadBase;
    QAction *actionSaveConfig;

    ArchArgs chipArgs;

    std::string currentBaseConfig;
};

NEXTPNR_NAMESPACE_END

#endif // MAINWINDOW_H
