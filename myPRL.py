###############################################################################
###############################################################################
##
##
##   myPRL vbeta - A simple Pressure Ruby Luminescence spectra fitting 
##   program for pressure determination in high pressure experiments.
##
##   Copyright (C) 2021 Alexis Forestier
##   E-mail : alforestier@gmail.com
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <https://www.gnu.org/licenses/>.
##
##
###############################################################################
###############################################################################

import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import pickle
import os
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from matplotlib.backend_bases import MouseButton
from scipy.spatial import ConvexHull
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


def lorentzian(x, *p):
    """ Lorentzian function, p is provided as (position, fwhm, area) 
        tuple argument with positional order."""
    x0, fwhm, A = p
    y = A*(2/(np.pi*fwhm))/(1 + ((x - x0)/(fwhm/2))**2)
    return y

def multi_lorentzian(x, *p):
    """ Sum of an arbitrary number of lorentzian, p provided as a longer tuple 
        of single lorentzian parameters, one triplet after the other. """
    y = 0
    for i in range(0, len(p), 3):
        _p = tuple(p[i:(i+3)])
        y += lorentzian(x, *_p)
    return y

def get_lorentzian_components(x, *p):
    """ Returns a DataFrame with each components of a multi_lorentzian sum 
    of lorentzian functions. p provided as a longer tuple 
        of single lorentzian parameters, one triplet after the other. """
    out = pd.DataFrame({'x':x})
    for i in range(0, len(p), 3):
        _p = tuple(p[i:(i+3)])
        y = pd.DataFrame({'y_'+str(int(i/3)):lorentzian(x, *_p)})
        out = out.join(y)
    return out


class Spec():
    """ A spectrum object containing raw data, 
    background and fits informations. """
    def __init__(self, filepath, header=[0], separator='\t', decimal='.',
                xunit = 'nm'):
        self.filepath = filepath
        self.filename = os.path.split(filepath)[1]
        self.data = pd.read_csv(self.filepath, header=header, 
                                sep=separator, decimal=decimal,
                                engine='python')
        self.data.columns = ['x','y']

        self.bg_points = None # background points
        self.fit_components = None
        self.fit_param = None

    def eval_rubber_band_bg(self):
        """ The rubberband background routine using ConvexHull to infer the 
        background of the spectrum"""
        x = np.array(self.data['x'])
        y = np.array(self.data['y'])
        v = ConvexHull(np.column_stack((x,y))).vertices
        v = np.roll(v, -v.argmin())
        anchors = v[:v.argmax()]
        self.bg_points = pd.DataFrame({'x':x[anchors], 'y':y[anchors]})
        self.data['bg'] = pd.Series(np.interp(x, x[anchors], y[anchors]))
        self.data['ycor'] = pd.Series(y - self.data['bg'])

    def eval_manual_bg(self):
        """ Manual background routine that fit the background on bg_points
        using spline"""
        x = np.array(self.data['x'])
        y = np.array(self.data['y'])

        # UnivariateSpline takes ascending ordered values:
        # bg_points must be set!
        x_bg = self.bg_points.sort_values(by='x')['x']
        y_bg = self.bg_points.sort_values(by='x')['y']

        # spline fit:
        spl_order = len(x_bg)-1 if len(x_bg)-1 <= 5 else 5
        spl = InterpolatedUnivariateSpline(x_bg, y_bg, k=spl_order)

        self.data['bg'] = pd.Series(spl(x))
        self.data['ycor'] = pd.Series(y - self.data['bg'])

    #def __del__(self):


class Ruby(Spec):
    """ A ruby luminescence spectrum with its own particular methods. """
    def __init__(self, lambda_laser=532, **kwargs):
        super().__init__(**kwargs)
        self.pressure = None
        self.lambda_ruby = None

    def eval_auto_fit(self, xunit='nm'):
        """ An autofit routine for two lorentzian components """
        x = np.array(self.data['x'])
        # use ycor unless no background was evaluated before
        if 'ycor' in self.data.columns:
            y = np.array(self.data['ycor'])
        else:
            y = np.array(self.data['y'])

        xbin = x[1] - x[0] # xunit/pt

        # min width of peaks 0.15 nm / 5 rel_cm_1 :
        pk, prop = find_peaks(y, height = max(y)/2, 
                    width = int(0.15/xbin) if xunit == 'nm' else int(5/xbin))

        param_ini = list()

        for i in range(len(pk)):
            # position:
            param_ini.append(x[pk[i]])
            # fwhm: 
            param_ini.append(x[int(prop['widths'][i])]-x[0])
            # area = pi*fwhm*I/2:
            param_ini.append(np.pi*(x[int(prop['widths'][i])]-x[0])*\
                                prop['peak_heights'][i]/2)

        # fit routine
        param, covar = curve_fit(multi_lorentzian, x, y, 
                                    p0=tuple(param_ini))

        # store fit results
        self.data['yfit'] = multi_lorentzian(x, *param)
        self.fit_components = get_lorentzian_components(x, *param)
        self.fit_param = param

    def eval_P(self, xunit='nm', lambda_laser=532, lambda_ruby_0=694.25):
        """ Pressure evaluation using the spectral 
        position of the 1st lorentzian peak. Ruby2020 scale is used."""
        if xunit == 'rel_cm_1':
            v_ruby_rel = self.fit_param[3]
            v_ruby = 10**7/lambda_laser - v_ruby_rel
            self.lambda_ruby = 10**7/v_ruby
        elif xunit == 'nm':
            self.lambda_ruby = self.fit_param[3]

        delta_lambda_ruby = self.lambda_ruby - lambda_ruby_0

        # ruby2020: High Pressure Research 
        # https://doi.org/10.1080/08957959.2020.1791107
        self.pressure = 1.87e3*(delta_lambda_ruby/lambda_ruby_0)* \
                (1 + 5.63*(delta_lambda_ruby/lambda_ruby_0))

    #def __del__(self):


class MyPRL(tk.Tk):
    """ Principal window of MyPRL with plot frame and action buttons"""
    def __init__(self):
        tk.Tk.__init__(self)
        self.geometry("1000x700")
        self.title("myPRL vbeta")
        self.resizable(True, True)
        self.columnconfigure(0, weight=2)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=10) # plot
        self.columnconfigure(3, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=10) # plot
        self.rowconfigure(1, weight=1)

        self._init_Vars()
        self._init_menubar()
        self._init_widgets()
        self._init_plot()

        self.rubies = {}
        self.rubygroups = {}

    def _init_Vars(self):
        #default options values:
        self.header = tk.IntVar(self, value='1')
        self.separator = tk.StringVar(self, value='\\t')
        self.decimal = tk.StringVar(self, value='.')
        self.xunit = tk.StringVar(self, value='nm')
        self.lambda_laser = tk.DoubleVar(self, value='532.00')
        self.lambda_ruby_0 = tk.DoubleVar(self, value='694.25')
        self.printedP = tk.DoubleVar(self, value='0.0')
        self.printedlambda = tk.DoubleVar(self, value='0.0')

        self._rubylist = tk.StringVar(self)

    def _init_menubar(self):
        menubar = tk.Menu(self)
        self.config(menu=menubar)

        menu_File = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=menu_File)
        menu_File.add_command(label="Load Files", command=self.load_files)
        menu_File.add_command(label='Load session...', 
                        command=self.load_session) 
        menu_File.add_command(label='Save session...', 
                        command=self.save_session) 
        menu_File.add_command(label="Quit", command=self.destroy)

    def _init_widgets(self):

        self.actionsframe = tk.LabelFrame(self, text="Actions")
        self.outputframe = tk.LabelFrame(self, text="Output")
        self.optionsframe = tk.LabelFrame(self, text="Global Options")

        # Left panel :
        self.button_load = tk.Button(self, text="Load Files", 
                                    command=self.load_files)
        self.button_remove = tk.Button(self.actionsframe, text="Remove", 
                                    command=self.remove_files)
        self.button_rubber_bg = tk.Button(self.actionsframe, 
                text="Rubberband Background", 
                command=lambda:[self.get_rubberband_bg(), self.plot_spec()],
                fg="blue")
        self.button_manual_bg = tk.Button(self.actionsframe, 
                                text="Manual Background", fg="blue",
                    command=lambda:[self.get_manual_bg(), self.plot_spec()])
        self.button_fit = tk.Button(self.actionsframe, text="Auto Fit",
                    command=lambda:[self.get_auto_fit_P(), self.plot_spec()],
                    fg="red")
        self.button_set_ref = tk.Button(self.actionsframe, 
                                text="Set as Ambient Ref.",
                                command=self.set_ruby_ref)
        self.button_group = tk.Button(self.actionsframe, 
                            text="Group",
                            command=lambda:[self.set_group(),self.plot_spec()], 
                            fg="green")

        self.filelist = tk.Listbox(self, height=25, 
                        listvariable=self._rubylist, selectmode = "extended",
                        activestyle="none")
        self.filelist.bind("<<ListboxSelect>>", self.plot_spec)


        self.scrollbar = tk.Scrollbar(self)
        self.scrollbar.config(command = self.filelist.yview)
        self.filelist.config(yscrollcommand = self.scrollbar.set)

        # Global options Frame :
        self.icon = ImageTk.PhotoImage(Image.open("./my_PRL_icon_large.png"))
        self.iconcanvas = tk.Canvas(self, width=220, height=180)
        self.iconcanvas.create_image(110, 90, 
                                    image=self.icon)

        self.label_fileformat = tk.Label(self.optionsframe, 
                                            text ="Files format:")
        self.label_header = tk.Label(self.optionsframe, 
                                        text ="Nb. of header line(s):")
        self.label_separator = tk.Label(self.optionsframe, text ="Separator:")
        self.label_decimal = tk.Label(self.optionsframe, text ="Decimal:")

        self.label_xunit = tk.Label(self.optionsframe, text ="x-axis unit:")
        self.label_lambda_laser = tk.Label(self.optionsframe, 
                                            text ="\u03BB Laser (nm):")

        self.entry_header = tk.Entry(self.optionsframe, 
                                        textvariable=self.header, width=8)
        self.entry_separator = tk.Entry(self.optionsframe, 
                                        textvariable=self.separator, width=8)
        self.entry_decimal = tk.Entry(self.optionsframe, 
                                        textvariable=self.decimal, width=8)

        self.entry_lambda_laser = tk.Entry(self.optionsframe, 
            textvariable=self.lambda_laser,
            state='disabled' if self.xunit.get() == 'nm' else 'normal', 
                                                            width=8)

        self.rdio_button_nm = tk.Radiobutton(self.optionsframe, 
            variable=self.xunit, 
            text='nm', value='nm',
            command = lambda: self.entry_lambda_laser.config(state='disabled'))
        self.rdio_button_relcm_1 = tk.Radiobutton(self.optionsframe, 
            variable=self.xunit, 
            text='rel. cm-1', value='rel_cm_1',
            command = lambda: self.entry_lambda_laser.config(state='normal'))

        self.label_ruby_scale = tk.Label(self.optionsframe, 
                                                        text ="Ruby scale:")
        self.label_ruby2020 = tk.Label(self.optionsframe, text ="Ruby2020")

        self.label_lambda_ruby_0 = tk.Label(self.optionsframe, 
                text ="\u03BB(P=0) (nm):")
        self.entry_lambda_ruby_0 = tk.Entry(self.optionsframe, 
                textvariable=self.lambda_ruby_0, width=8)

        self.separatorline1 = ttk.Separator(self.optionsframe, 
                                    orient='horizontal')
        self.separatorline2 = ttk.Separator(self.optionsframe, 
                                    orient='horizontal')


        # output widgets:
        self.labelprintedP = tk.Label(self.outputframe, 
                text ="P (GPa) = ")
        self.labelprintedlambda = tk.Label(self.outputframe, 
                text ="\u03BB (nm) = ")
        self.output_printedP = tk.Entry(self.outputframe, 
                textvariable=self.printedP, width=10)
        self.output_printedlambda = tk.Entry(self.outputframe, 
                textvariable=self.printedlambda, width=10)
        self.button_export_P = tk.Button(self.outputframe, 
                            text="Export Pressure Values",
                            command=self.export_P)

        # grid/pack :

        self.optionsframe.grid(row=1, column=3, sticky='news')
        self.iconcanvas.grid(row=2, column=3, sticky='news')

        # Global Options Frame:
        self.label_fileformat.grid(row=0, column=0, sticky='w')

        self.label_header.grid(row=1, column=0, sticky='e')
        self.entry_header.grid(row=1, column=1, sticky='ew')

        self.label_separator.grid(row=2, column=0, sticky='e')
        self.entry_separator.grid(row=2, column=1, sticky='ew')

        self.label_decimal.grid(row=3, column=0, sticky='e')
        self.entry_decimal.grid(row=3, column=1, sticky='ew')

        self.separatorline1.grid(row=4, columnspan=2, sticky='ew')
        self.label_xunit.grid(row=5, column=0, sticky='w')
        self.rdio_button_nm.grid(row=6, column=0, sticky='news')
        self.rdio_button_relcm_1.grid(row=6, column=1, sticky='news')

        self.label_lambda_laser.grid(row=7, column=0, sticky='e')
        self.entry_lambda_laser.grid(row=7, column=1, sticky='ew')

        self.separatorline2.grid(row=8, columnspan=2, sticky='ew')

        self.label_ruby_scale.grid(row=9, column=0, sticky='e')
        self.label_ruby2020.grid(row=9, column=1, sticky='ew')

        self.label_lambda_ruby_0.grid(row=10, column=0, sticky='e')
        self.entry_lambda_ruby_0.grid(row=10, column=1, sticky='ew')

        # actions frame
        self.actionsframe.grid(row=2, column=0, sticky='news')

        self.button_load.grid(row=0, column=0, sticky='news')
        self.filelist.grid(row=1, column=0, sticky='news')
        self.scrollbar.grid(row=1, column=1, sticky='nsw')
        self.button_rubber_bg.pack(fill='both')
        self.button_manual_bg.pack(fill='both')
        self.button_fit.pack(fill='both')
        self.button_group.pack(fill='both')
        self.button_set_ref.pack(fill='both')
        self.button_remove.pack(fill='both')

        # output frame:
        self.outputframe.grid(row=2, column=2, sticky='news')
        self.labelprintedP.grid(row=0, column=0, sticky='ne')
        self.output_printedP.grid(row=0, column=1, sticky='ne')
        self.labelprintedlambda.grid(row=1, column=0, sticky='ne')
        self.output_printedlambda.grid(row=1, column=1, sticky='ne')
        self.button_export_P.grid(row=2, columnspan=2, sticky='ne')

    def _init_plot(self):
        self.fig = Figure()
        self.plot = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)

        self.toolbarFrame = tk.Frame()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        #self.toolbar = NavigationToolbar2Tk(self.canvas, self, 
        #                                    pack_toolbar=False)
        #self.toolbar.grid(row=0, column=2, sticky='news')
        self.toolbarFrame.grid(row=0, column=2, sticky='news')
        self.canvas.get_tk_widget().grid(row=1, column=2, sticky='news')

        # will be used as ax matplotlib object for on the fly background in 
        # manual mode:
        self._otf_bg = None

    def _get_selected_spec(self):
        """ Returns filename list corresponding 
        to selection in filelist."""
        selected = [self.filelist.get(elem) \
                    for elem in self.filelist.curselection()]

        return selected

    def _get_spec_listboxpos(self, names):
        """ Returns a list of integer corresponding to 
            positions in the filelist listbox 
            corresponding to the names input."""
        pos = [i for i in range(self.filelist.size()) \
                    if self.filelist.get(i) in names]
        return pos

    def _update_rubylist_Var(self):
        """ Automatic update of the ListBox StringVar corresponding 
        to the current rubies dict. state."""
        self._rubylist.set(list(self.rubies.keys()))

    def _update_printed_P_lambda(self, P, lambda_ruby):
        """ update printed values in output. """
        self.printedP.set(P)
        self.printedlambda.set(lambda_ruby)

    def load_files(self):
        """ Load and store multiple files in rubies object"""
        filepaths = tk.filedialog.askopenfilenames(title="Load Files",
                                        filetypes=[('all files','.*')])
        try:
            for path in filepaths:
                ruby = Ruby(filepath=path, 
                  header=list(range(self.header.get())) \
                    if self.header.get() != 0 else None,
                    separator=str(self.separator.get()),
                    decimal=str(self.decimal.get()),
                    xunit=str(self.xunit.get()),
                    lambda_laser=self.lambda_laser.get())            

                if ruby.filename not in self.rubies:
                    # fill rubies dict:
                    self.rubies[ruby.filename] = ruby 
                else:
                    tk.messagebox.showwarning("Warning !",
                                    "File {} already exists.\n"\
                                    .format(ruby.filename) +
                                    "It has not been replaced !")

            self._update_rubylist_Var()

        except (ValueError, pd.errors.ParserError) as e:
            tk.messagebox.showerror("Loading failed", 
                                   "Loading failed: \ncheck data file format "
                                   "in File>Options\n\n"
                                   "Details: " + str(e))

    def remove_files(self):
        """ Remove selected files from the software"""
        selec = self._get_selected_spec()

        for name in selec:
            del self.rubies[name]
            # del groups:
            if name in self.rubygroups: 
                del self.rubygroups[name]

        self._update_rubylist_Var()

    def get_rubberband_bg(self):
        """ Get background for all spectra selected """
        selec = self._get_selected_spec()

        try:
            for name in selec:
                self.rubies[name].eval_rubber_band_bg()
        except Exception as e:
            tk.messagebox.showerror("Rubberband background routine failed", 
                                      "Rubberband background routine "
                                      "failed...\n\n"
                                      "Details: " + str(e))

    def get_manual_bg(self):
        """ Manual background mode, one file only ! """

        selec = self._get_selected_spec()

        # activate manual background mode:
        if self.button_manual_bg.cget('bg') != 'blue':

            if(len(selec) > 1):
                tk.messagebox.showwarning("Warning",
                    "Manual background routine run over one data file only.\n"
                    "Only the first selection will be affected")

            # only keep the first selection:
            selec = selec[0]

            try:
                # reset background:
                self.rubies[selec].bg_points = pd.DataFrame()
                if 'bg' in self.rubies[selec].data:
                    del self.rubies[selec].data['bg']
                if 'ycor' in self.rubies[selec].data:
                    del self.rubies[selec].data['ycor']

                self.button_manual_bg.configure(bg='blue', 
                                            activebackground='blue',
                                            fg="white")
                self.fig.suptitle("Manual Background Mode > Use Right Click", 
                                size='xx-large',
                                color='blue', weight="bold")

                self.fig.canvas.mpl_connect('button_press_event', 
                                        self.set_manual_bg_points)

            except Exception as e: 
                tk.messagebox.showerror("Manual background routine failed", 
                                    "Manual background routine failed...\n\n"
                                    "Details: " + str(e))

        # deactivate manual background mode:
        else:
            self.button_manual_bg.configure(bg='#d9d9d9', 
                                            activebackground='#d9d9d9',
                                            fg="blue")
            self.fig.suptitle("")
            self.fig.canvas.mpl_disconnect(
                            self.fig.canvas.mpl_connect('button_press_event', 
                            self.set_manual_bg_points))
            self.plot.autoscale(True)

    def set_manual_bg_points(self, event):
        """ Routine to add one individual point
            in bg_points DataFrame attribute of Spec object"""
 
        # only keep the first selection:
        selec = self._get_selected_spec()[0]

        # keep autoscale disabled (draw is called each time!)
        self.plot.autoscale(False)

        if event.button is MouseButton.RIGHT:

            point = pd.DataFrame({'x':[event.xdata], 'y':[event.ydata]})

            self.rubies[selec].bg_points = \
            self.rubies[selec].bg_points.append(point, ignore_index=True)

            self.plot.scatter(point['x'], point['y'], 
                    color="darkblue", s=250, marker="+")

            x = self.rubies[selec].data['x']

            if len(self.rubies[selec].bg_points) == 2:
                self.rubies[selec].eval_manual_bg()

                bg = self.rubies[selec].data['bg']
                # create the on the fly background artist:
                self._otf_bg, = self.plot.plot(x, bg, 'b--')

            elif len(self.rubies[selec].bg_points) > 2:
                self.rubies[selec].eval_manual_bg()

                bg = self.rubies[selec].data['bg']
                # update the on the fly background artist:
                self._otf_bg.set_data(x, bg)

            self.canvas.draw()

    def get_auto_fit_P(self):
        """ Perform autofit for all spectra selected """
        selec = self._get_selected_spec()

        try:
            for name in selec:
                self.rubies[name].eval_auto_fit(xunit=str(self.xunit.get()))
                self.rubies[name].eval_P(xunit=str(self.xunit.get()),
                                    lambda_laser=self.lambda_laser.get(),
                                        lambda_ruby_0=self.lambda_ruby_0.get())
        except Exception as e:
            tk.messagebox.showerror("Autofit failed", 
                                      "Autofit failed...\n\n"
                                      "Details: " + str(e))

    def set_ruby_ref(self):
        """ Use the current ruby spectrum selection as ambient pressure 
        reference. """
        try:
            # keep only the first file if multiple selection is used:
            selec = self._get_selected_spec()[0]

            l_ruby_0_list = [self.rubies[selec].lambda_ruby]
            # use average wavelength if grouped:
            if selec in self.rubygroups:
                for friend in self.rubygroups[selec]:
                    l_ruby_f = self.rubies[friend].lambda_ruby
                    l_ruby_0_list.append(l_ruby_f)

            l_ruby_0 = sum(l_ruby_0_list)/len(l_ruby_0_list)
            self.lambda_ruby_0.set(round(l_ruby_0, 4))

            # update pressure values already fitted:
            for name in self.rubies.keys():
                if self.rubies[name].pressure != None:
                    self.rubies[name].eval_P(xunit=str(self.xunit.get()),
                                    lambda_laser=self.lambda_laser.get(),
                                    lambda_ruby_0=self.lambda_ruby_0.get())
            # replot:
            self.plot_spec()

        except Exception as e:
            tk.messagebox.showerror("Set ruby ref. failed", 
                                      "Set ruby ref. failed...\n\n"
                                      "Details: " + str(e))
    def set_group(self):
        """ Group ruby spectra in order to use average pressure """
        selec = self._get_selected_spec()

        # create group
        if self.button_group.cget('bg') == '#d9d9d9':
            if(len(selec) >= 2):
                for name in selec:
                    # all selected files except itself:
                    friends = [x for x in selec if x != name]
                    self.rubygroups[name] = friends

                    self.button_group.configure(bg='green', 
                                            activebackground='green',
                                            fg="white",
                                            text="Ungroup")
            else:
                pass

        # delete group
        elif self.button_group.cget('bg') == 'green':
            for name in selec:
                if name in self.rubygroups: 
                    del self.rubygroups[name]
            self.filelist.selection_clear(0, 'end')
            self.filelist.selection_set(self._get_spec_listboxpos(selec[0]))

            self.button_group.configure(bg='#d9d9d9', 
                                            activebackground='#d9d9d9',
                                            fg="black",
                                            text="Group")

    def _group_selec(self, selec):
        """ Check if grouped or not and change the state of the group button 
        accordingly. 
        Set selection of the whole group in filelist, and return the list
        of rubies grouped to add to the plot. 
        """
        if list(set(self.rubygroups.keys()).intersection(selec)) != []:
            self.button_group.configure(bg='green', 
                                        activebackground='green',
                                        fg="white",
                                        text="Ungroup")
        else:
            self.button_group.configure(bg='#d9d9d9', 
                                        activebackground='#d9d9d9',
                                        fg="black",
                                        text="Group")
        sel = []
        for name in selec:
            if name in self.rubygroups: 
                friends = self.rubygroups[name]
            else:
                friends = []
            for f in friends:
                if f not in selec:
                    pos = self._get_spec_listboxpos(f)
                    self.filelist.selection_set(pos)
                    sel.append(f)
        return sel

    def plot_spec(self, *args):
        """ The plot spectra routine. """
        selec = self._get_selected_spec()
        selec += self._group_selec(selec)

        self.plot.cla() # clear for the new plot

        for name in selec:
            x = self.rubies[name].data['x']
            y = self.rubies[name].data['y']
            self.plot.scatter(x,y,color="gray", alpha=.1, s=20)  

            # if background exist then plot it !
            if 'bg' in self.rubies[name].data:
                bg_points = self.rubies[name].bg_points # is a pd.DataFrame
                bg = self.rubies[name].data['bg']
                ycor = self.rubies[name].data['ycor']
                self.plot.scatter(bg_points['x'], bg_points['y'], 
                    color="darkblue", s=250, marker="+")
                self.plot.plot(x,bg,'b--')
                self.plot.scatter(x, ycor, color="black", alpha=.6, s=20)

            # if fit exist then plot it !
            if 'yfit' in self.rubies[name].data:
                self.plot.plot(x, self.rubies[name].data['yfit'], 
                    '-', color='darkorange', label='best fit')
                self.plot.plot(x, self.rubies[name].fit_components.iloc[:,1:], 
                    '--', color='darkred', label='component')
                self.plot.axvline(x = self.rubies[name].fit_param[3],
                     color='k')

                text_x = self.rubies[name].fit_param[3]
                # max of the ruby line:
                text_y = 0.9*(2*self.rubies[name].fit_param[5]/(np.pi*\
                        self.rubies[name].fit_param[4]))

                # print pressure value in the plot:
                if name not in self.rubygroups:
                    p = self.rubies[name].pressure
                    l = self.rubies[name].lambda_ruby

                    self.plot.text(x=text_x, 
                                   y=text_y, 
                                   s= '    \u03BB = ' + 
                                   str(round(l,4)) + ' nm' +
                                   '\n    P = ' +
                                   str(round(p,2)) + ' GPa',
                                   size=12, weight='bold')

                    self._update_printed_P_lambda(round(p,2),round(l,4))

                else:
                    plist = [self.rubies[name].pressure] + \
                    [self.rubies[k].pressure for k in self.rubygroups[name]]
                    # average:
                    p = sum(plist)/len(plist)

                    self.plot.text(x=text_x, 
                                   y=text_y, 
                                   s= '   Pm = ' + \
                                   str(round(p,2)) + ' GPa',
                                   size=12, weight='bold',
                                   color='green')  

                    self._update_printed_P_lambda(round(p,2),'NA')

        self.canvas.draw()
        self.toolbar.update()

    def export_P(self):
        """ Export calculated pressure values in a txt file. """
        
        path = tk.filedialog.asksaveasfilename(
                    filetypes=[('Text file', '.txt')], 
                    defaultextension='.txt')

        values = pd.DataFrame({'file':(),'pressure':()})

        for name in self.rubies:

            already = [name for group in list(values['file']) \
             for name in group]

            if self.rubies[name].pressure != None and name not in already:
                _P = {'file':(name,), 'pressure':(self.rubies[name].pressure,)}
                if name in self.rubygroups:
                    for friend in self.rubygroups[name]:
                        _P['file'] += (friend,)
                        _P['pressure'] += (self.rubies[friend].pressure,)
            else:
                _P = None

            if _P != None:
                values = values.append(_P, ignore_index=True)

        export = pd.DataFrame({'file':
                    list(map(lambda s: ', '.join(s), values["file"])),
                    'pressure':
                    list(map(lambda s: sum(s)/len(s), values["pressure"]))})

        # file = open(path, 'w')
        # file.write(export.to_string())
        # file.close()
        export.to_csv(path, sep='\t', header=True, index=False)

    def save_session(self):
        """ Save current session in a binary *.pkl 
        file that can be loaded later. """

        file_path = tk.filedialog.asksaveasfilename(
                    filetypes=[('Binary file (pickle)', '.pkl')], 
                    defaultextension='.pkl')

        options = {'header':self.header.get(),
                   'separator':self.separator.get(),
                   'decimal':self.decimal.get(),
                   'xunit':self.xunit.get(),
                   'lambda_laser':self.lambda_laser.get(),
                   'lambda_ruby_0':self.lambda_ruby_0.get()}

        try:
            file = open(file_path,'wb')
            pickle.dump(options, file)
            pickle.dump(self.rubies, file)
            pickle.dump(self.rubygroups, file)
            file.close()
        except Exception as e:
            tk.messagebox.showerror("Save session routine failed", 
                                    "Save session routine failed...\n\n"
                                    "Details: " + str(e))
    def load_session(self):

        file_path = tk.filedialog.askopenfilename(
                    filetypes=[('Binary file (pickle)', '.pkl')], 
                    defaultextension='.pkl')

        try:
            file = open(file_path, 'rb')
            loaded_options = pickle.load(file)
            loaded_rubies = pickle.load(file)
            loaded_rubygroups = pickle.load(file)
            file.close()
    
            # load options:
            self.header.set(loaded_options['header'])
            self.separator.set(loaded_options['separator'])
            self.decimal.set(loaded_options['decimal'])
            self.xunit.set(loaded_options['xunit'])
            self.lambda_laser.set(loaded_options['lambda_laser'])
            self.lambda_ruby_0.set(loaded_options['lambda_ruby_0'])
                
            # load rubies:
            for name in loaded_rubies:
                if name not in self.rubies:
                    self.rubies[name] = loaded_rubies[name]
                else:
                    tk.messagebox.showerror("Load session routine failed",
                                    "File {} already exists.\n"\
                                    .format(name) +
                                    "It has not been replaced !")

            self._update_rubylist_Var()

            #load rubygroups:
            self.rubygroups = loaded_rubygroups

        except Exception as e:
            tk.messagebox.showerror("Load session routine failed", 
                                     "Load session routine failed...\n\n"
                                     "Details: " + str(e))

if __name__ == '__main__':
    root = MyPRL()
    root.mainloop()

