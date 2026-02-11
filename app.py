import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class CTViewerGUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("CT Viewer GUI")
        self.geometry("1100x750")
        self.minsize(900, 650)

        # value: dict {path, data, shape, meta}
        self.volumes = {}
        self.current_key = None

        # Slider ranges default
        self.axial_min, self.axial_max = 0, 10
        self.sagital_min, self.sagital_max = 0, 10
        self.coronal_min, self.coronal_max = 0, 10

        self.axial_var = tk.IntVar(value=0)
        self.sagital_var = tk.IntVar(value=0)
        self.coronal_var = tk.IntVar(value=0)

        self._build_styles()
        self._build_menu_bar()
        self._build_top_bar()
        self._build_main_layout()

        # Mapping view -> (axes, canvas) untuk kontrol zoom
        self.view_map = {
            "axial": (self.axial_axes, self.axial_canvas),
            "sagital": (self.sagital_axes, self.sagital_canvas),
            "coronal": (self.coronal_axes, self.coronal_canvas),
        }

        # Simpan zoom state per view
        self.zoom_state = {
            "axial": {"base": None, "current": None},
            "sagital": {"base": None, "current": None},
            "coronal": {"base": None, "current": None},
        }

        self._render_placeholders()
        self._set_metadata_rows({})  # default tampilkan placeholder metadata
        self._show_left_mode("Metadata")  # default mode

    def _build_styles(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure("TFrame", background="white")
        style.configure("App.TFrame", background="white")
        style.configure("Top.TFrame", background="white")
        style.configure("Panel.TLabelframe", background="white")
        style.configure("Panel.TLabelframe.Label", background="white", font=("Segoe UI", 10, "bold"))
        style.configure("TLabel", background="white", font=("Segoe UI", 10))
        style.configure("TButton", font=("Segoe UI", 10))
        style.configure("TCombobox", font=("Segoe UI", 10))
        style.configure("Treeview", font=("Segoe UI", 9))
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))

        self.configure(bg="white")

    # =========================
    # MENU BAR
    # =========================
    def _build_menu_bar(self):
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Import NIfTI (.nii atau .nii.gz)", command=self.import_nifti)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.destroy)

        menubar.add_cascade(label="File", menu=file_menu)
        self.config(menu=menubar)

    # =========================
    # TOP BAR
    # =========================
    def _build_top_bar(self):
        top = ttk.Frame(self, style="Top.TFrame", padding=(16, 16, 16, 8))
        top.grid(row=0, column=0, sticky="ew")
        self.grid_columnconfigure(0, weight=1)

        ttk.Label(top, text="DROPDOWN DATASET").grid(row=0, column=0, sticky="w", padx=(0, 8))

        self.dataset_var = tk.StringVar(value="Belum ada file")
        self.dataset_dropdown = ttk.Combobox(
            top,
            textvariable=self.dataset_var,
            values=["Belum ada file"],
            state="readonly",
            width=40
        )
        self.dataset_dropdown.grid(row=0, column=1, sticky="w")
        self.dataset_dropdown.bind("<<ComboboxSelected>>", self.on_select_dataset)

        top.grid_columnconfigure(2, weight=1)

        self.save_btn = ttk.Button(top, text="TOMBOL SAVE", command=self.on_save)
        self.save_btn.grid(row=0, column=3, sticky="e")

    # =========================
    # MAIN LAYOUT
    # =========================
    def _build_main_layout(self):
        container = ttk.Frame(self, style="App.TFrame", padding=(16, 8, 16, 16))
        container.grid(row=1, column=0, sticky="nsew")
        self.grid_rowconfigure(1, weight=1)

        container.grid_columnconfigure(0, weight=0)   # kiri
        container.grid_columnconfigure(1, weight=1)   # kanan
        container.grid_rowconfigure(0, weight=1)

        # ===========================
        # KIRI: PANEL DENGAN DROPDOWN MODE
        # ===========================
        left_panel = ttk.Labelframe(
            container,
            text="PANEL KIRI",
            style="Panel.TLabelframe",
            padding=(10, 10)
        )
        left_panel.grid(row=0, column=0, sticky="nsew", padx=(0, 12))
        left_panel.configure(width=320)
        left_panel.grid_propagate(False)

        left_panel.grid_rowconfigure(1, weight=1)
        left_panel.grid_columnconfigure(0, weight=1)

        # Dropdown mode (Metadata / Segmentasi)
        mode_bar = ttk.Frame(left_panel, style="App.TFrame")
        mode_bar.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        mode_bar.grid_columnconfigure(1, weight=1)

        ttk.Label(mode_bar, text="Mode").grid(row=0, column=0, sticky="w", padx=(0, 8))

        self.left_mode_var = tk.StringVar(value="Metadata")
        self.left_mode_dropdown = ttk.Combobox(
            mode_bar,
            textvariable=self.left_mode_var,
            values=["Metadata", "Segmentasi"],
            state="readonly",
            width=14
        )
        self.left_mode_dropdown.grid(row=0, column=1, sticky="w")
        self.left_mode_dropdown.bind("<<ComboboxSelected>>", self.on_left_mode_change)

        # Container untuk menampung 2 halaman (metadata dan segmentasi)
        self.left_pages = ttk.Frame(left_panel, style="App.TFrame")
        self.left_pages.grid(row=1, column=0, sticky="nsew")
        self.left_pages.grid_rowconfigure(0, weight=1)
        self.left_pages.grid_columnconfigure(0, weight=1)

        # ==== Page: Metadata ====
        self.meta_page = ttk.Frame(self.left_pages, style="App.TFrame")
        self.meta_page.grid(row=0, column=0, sticky="nsew")
        self.meta_page.grid_rowconfigure(1, weight=1)
        self.meta_page.grid_columnconfigure(0, weight=1)

        ttk.Label(
            self.meta_page,
            text="Informasi dari header dan properti volume. Jika tidak ada nilai, ditampilkan sebagai -.",
            wraplength=280,
            justify="left"
        ).grid(row=0, column=0, sticky="ew", pady=(0, 10))

        self.meta_tree = ttk.Treeview(
            self.meta_page,
            columns=("field", "value"),
            show="headings",
            height=20
        )
        self.meta_tree.heading("field", text="Field")
        self.meta_tree.heading("value", text="Value")
        self.meta_tree.column("field", width=140, anchor="w")
        self.meta_tree.column("value", width=160, anchor="w")
        self.meta_tree.grid(row=1, column=0, sticky="nsew")

        meta_scroll = ttk.Scrollbar(self.meta_page, orient="vertical", command=self.meta_tree.yview)
        self.meta_tree.configure(yscrollcommand=meta_scroll.set)
        meta_scroll.grid(row=1, column=1, sticky="ns")

        # ==== Page: Segmentasi (placeholder) ====
        self.seg_page = ttk.Frame(self.left_pages, style="App.TFrame")
        self.seg_page.grid(row=0, column=0, sticky="nsew")
        self.seg_page.grid_rowconfigure(1, weight=1)
        self.seg_page.grid_columnconfigure(0, weight=1)

        ttk.Label(
            self.seg_page,
            text="Tool Segmentasi (placeholder)\n\nSilakan isi:\n• Brush\n• Eraser\n• Threshold\n• Undo/Redo\n• Overlay Mask\n• Export Mask",
            justify="left",
            wraplength=280
        ).grid(row=0, column=0, sticky="nw")

        # Default: metadata page tampil
        self.meta_page.tkraise()

        # ===========================
        # KANAN: VISUALISASI
        # ===========================
        vis = ttk.Frame(container, style="App.TFrame")
        vis.grid(row=0, column=1, sticky="nsew")

        vis.grid_columnconfigure(0, weight=1)
        vis.grid_columnconfigure(1, weight=1)
        vis.grid_rowconfigure(0, weight=0)
        vis.grid_rowconfigure(1, weight=1)
        vis.grid_rowconfigure(2, weight=0)
        vis.grid_rowconfigure(3, weight=1)

        # sliders axial dan sagital
        self.axial_slider_row = self._build_slider_row(
            parent=vis,
            label_text="SLIDER UNTUK GANTI SLICE AXIAL",
            var=self.axial_var,
            store_scale_as="axial_scale",
            on_change=self.on_axial_change
        )
        self.axial_slider_row.grid(row=0, column=0, sticky="ew", padx=4, pady=(0, 6))

        self.sagital_slider_row = self._build_slider_row(
            parent=vis,
            label_text="SLIDER UNTUK GANTI SLICE SAGITAL",
            var=self.sagital_var,
            store_scale_as="sagital_scale",
            on_change=self.on_sagital_change
        )
        self.sagital_slider_row.grid(row=0, column=1, sticky="ew", padx=4, pady=(0, 6))

        # views axial dan sagital
        self.axial_panel, self.axial_axes, self.axial_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI AXIAL",
            view_key="axial"
        )
        self.axial_panel.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        self.sagital_panel, self.sagital_axes, self.sagital_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI SAGITAL",
            view_key="sagital"
        )
        self.sagital_panel.grid(row=1, column=1, sticky="nsew", padx=4, pady=4)

        # slider coronal
        self.coronal_slider_row = self._build_slider_row(
            parent=vis,
            label_text="SLIDER UNTUK GANTI SLICE CORONAL",
            var=self.coronal_var,
            store_scale_as="coronal_scale",
            on_change=self.on_coronal_change
        )
        self.coronal_slider_row.grid(row=2, column=0, sticky="ew", padx=4, pady=(6, 6))
        ttk.Frame(vis, style="App.TFrame").grid(row=2, column=1, sticky="ew")

        # view coronal dan 3D
        self.coronal_panel, self.coronal_axes, self.coronal_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI CORONAL",
            view_key="coronal"
        )
        self.coronal_panel.grid(row=3, column=0, sticky="nsew", padx=4, pady=4)

        panel3d = ttk.Labelframe(vis, text="GRID UNTUK 3D", style="Panel.TLabelframe", padding=(10, 10))
        panel3d.grid(row=3, column=1, sticky="nsew", padx=4, pady=4)
        panel3d.grid_rowconfigure(0, weight=1)
        panel3d.grid_columnconfigure(0, weight=1)
        ttk.Label(panel3d, text="TAMPILAN 3D (placeholder)", anchor="center").grid(row=0, column=0, sticky="nsew")

    # =========================
    # LEFT MODE SWITCH
    # =========================
    def on_left_mode_change(self, event=None):
        self._show_left_mode(self.left_mode_var.get())

    def _show_left_mode(self, mode: str):
        if mode == "Segmentasi":
            self.seg_page.tkraise()
        else:
            self.meta_page.tkraise()

    # =========================
    # SLIDER ROW
    # =========================
    def _build_slider_row(self, parent, label_text, var, store_scale_as, on_change):
        frame = ttk.Frame(parent, style="App.TFrame")
        frame.grid_columnconfigure(2, weight=1)

        ttk.Label(frame, text=label_text).grid(row=0, column=0, sticky="w", padx=(0, 10))

        btn_left = ttk.Button(frame, text="◀", width=3, command=lambda: self._step(var, -1, on_change))
        btn_left.grid(row=0, column=1, sticky="w")

        slider = ttk.Scale(
            frame,
            from_=0,
            to=10,
            orient="horizontal",
            command=lambda val: self._set_from_slider(var, val, on_change)
        )
        slider.grid(row=0, column=2, sticky="ew", padx=(8, 8))

        btn_right = ttk.Button(frame, text="▶", width=3, command=lambda: self._step(var, 1, on_change))
        btn_right.grid(row=0, column=3, sticky="e")

        value_lbl = ttk.Label(frame, text=str(var.get()), width=5, anchor="e")
        value_lbl.grid(row=0, column=4, sticky="e", padx=(10, 0))

        setattr(self, store_scale_as, slider)

        def _sync(*_):
            v = var.get()
            min_v = int(float(slider.cget("from")))
            max_v = int(float(slider.cget("to")))
            if v < min_v:
                v = min_v
                var.set(v)
            if v > max_v:
                v = max_v
                var.set(v)
            value_lbl.config(text=str(v))
            slider.set(v)

        var.trace_add("write", _sync)
        return frame

    # =========================
    # VIEW ONLY + ZOOM BUTTONS
    # =========================
    def _build_view_plot_only(self, parent, title, view_key):
        panel = ttk.Labelframe(parent, text=title, style="Panel.TLabelframe", padding=(6, 6))
        panel.grid_rowconfigure(1, weight=1)
        panel.grid_columnconfigure(0, weight=1)

        tools = ttk.Frame(panel, style="App.TFrame")
        tools.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        tools.grid_columnconfigure(0, weight=1)

        btn_fit = ttk.Button(tools, text="Fit", width=6, command=lambda: self.zoom_fit(view_key))
        btn_out = ttk.Button(tools, text="Zoom Out", width=9, command=lambda: self.zoom_out(view_key))
        btn_in = ttk.Button(tools, text="Zoom In", width=8, command=lambda: self.zoom_in(view_key))

        btn_in.pack(side="right", padx=(6, 0))
        btn_out.pack(side="right", padx=(6, 0))
        btn_fit.pack(side="right")

        fig = Figure(figsize=(4, 4), dpi=100)
        ax = fig.add_subplot(111)
        ax.set_axis_off()

        canvas = FigureCanvasTkAgg(fig, master=panel)
        canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")

        return panel, ax, canvas

    # =========================
    # METADATA HELPERS
    # =========================
    @staticmethod
    def _fmt(value):
        if value is None:
            return "-"
        if isinstance(value, (list, tuple)):
            if len(value) == 0:
                return "-"
            return str(value)
        if isinstance(value, np.ndarray):
            if value.size == 0:
                return "-"
            return np.array2string(value, precision=4, separator=", ")
        s = str(value).strip()
        return s if s else "-"

    def _set_metadata_rows(self, meta: dict):
        for item in self.meta_tree.get_children():
            self.meta_tree.delete(item)

        if not meta:
            default_fields = [
                "FileName", "FilePath", "Format",
                "Dim", "Shape", "VoxelSize",
                "DataType", "BitPix",
                "SFormCode", "QFormCode",
                "Affine",
                "CalMin", "CalMax",
                "Descrip"
            ]
            for f in default_fields:
                self.meta_tree.insert("", "end", values=(f, "-"))
            return

        for k in sorted(meta.keys(), key=lambda x: x.lower()):
            self.meta_tree.insert("", "end", values=(k, self._fmt(meta.get(k))))

    def _extract_nifti_metadata(self, img, data, path: str) -> dict:
        meta = {}
        meta["FileName"] = os.path.basename(path)
        meta["FilePath"] = os.path.abspath(path)
        meta["Format"] = "NIfTI"
        meta["Dim"] = int(data.ndim)
        meta["Shape"] = tuple(int(i) for i in data.shape)

        hdr = img.header
        try:
            meta["DataType"] = str(hdr.get_data_dtype())
        except Exception:
            meta["DataType"] = None

        try:
            meta["BitPix"] = int(hdr["bitpix"])
        except Exception:
            meta["BitPix"] = None

        try:
            zooms = hdr.get_zooms()
            if zooms is not None and len(zooms) >= 3:
                meta["VoxelSize"] = (float(zooms[0]), float(zooms[1]), float(zooms[2]))
            else:
                meta["VoxelSize"] = None
        except Exception:
            meta["VoxelSize"] = None

        try:
            des = hdr.get("descrip", b"")
            if isinstance(des, (bytes, bytearray)):
                des = des.decode(errors="ignore")
            meta["Descrip"] = des
        except Exception:
            meta["Descrip"] = None

        try:
            meta["CalMin"] = float(hdr.get("cal_min", np.nan))
            if np.isnan(meta["CalMin"]):
                meta["CalMin"] = None
        except Exception:
            meta["CalMin"] = None

        try:
            meta["CalMax"] = float(hdr.get("cal_max", np.nan))
            if np.isnan(meta["CalMax"]):
                meta["CalMax"] = None
        except Exception:
            meta["CalMax"] = None

        try:
            meta["QFormCode"] = int(hdr["qform_code"])
        except Exception:
            meta["QFormCode"] = None

        try:
            meta["SFormCode"] = int(hdr["sform_code"])
        except Exception:
            meta["SFormCode"] = None

        try:
            meta["Affine"] = np.asarray(img.affine, dtype=np.float64)
        except Exception:
            meta["Affine"] = None

        return meta

    # =========================
    # IMPORT DAN SELEKSI DATASET
    # =========================
    def import_nifti(self):
        path = filedialog.askopenfilename(
            title="Pilih file NIfTI",
            filetypes=[("NIfTI files", "*.nii *.nii.gz"), ("All files", "*.*")]
        )
        if not path:
            return

        try:
            import nibabel as nib
        except ModuleNotFoundError:
            messagebox.showerror(
                "Nibabel belum ada",
                "Install nibabel di environment yang sama:\n\npython -m pip install nibabel"
            )
            return

        try:
            img = nib.load(path)
            data = img.get_fdata()
        except Exception as e:
            messagebox.showerror("Gagal load NIfTI", f"Error:\n{e}")
            return

        if data.ndim == 4:
            data = data[..., 0]

        if data.ndim != 3:
            messagebox.showerror("Dimensi tidak didukung", f"Shape data: {data.shape}")
            return

        data = np.asarray(data, dtype=np.float32)
        meta = self._extract_nifti_metadata(img, data, path)

        base = os.path.basename(path)
        key = self._make_unique_key(base)

        self.volumes[key] = {"path": path, "data": data, "shape": data.shape, "meta": meta}

        self._refresh_dataset_dropdown(select_key=key)
        self._set_current_dataset(key)

    def _make_unique_key(self, base_name: str) -> str:
        if base_name not in self.volumes:
            return base_name
        i = 2
        while True:
            candidate = f"{base_name} ({i})"
            if candidate not in self.volumes:
                return candidate
            i += 1

    def _refresh_dataset_dropdown(self, select_key=None):
        keys = list(self.volumes.keys())
        if not keys:
            keys = ["Belum ada file"]

        self.dataset_dropdown.configure(values=keys)

        if select_key is not None and select_key in self.volumes:
            self.dataset_var.set(select_key)
        elif self.current_key in self.volumes:
            self.dataset_var.set(self.current_key)
        else:
            self.dataset_var.set(keys[0])

    def on_select_dataset(self, event=None):
        key = self.dataset_var.get()
        if key in self.volumes:
            self._set_current_dataset(key)

    def _set_current_dataset(self, key: str):
        self.current_key = key
        vol = self.volumes[key]
        data = vol["data"]

        # Update metadata panel
        self._set_metadata_rows(vol.get("meta", {}))

        x, y, z = data.shape

        self.axial_min, self.axial_max = 0, z - 1
        self.sagital_min, self.sagital_max = 0, x - 1
        self.coronal_min, self.coronal_max = 0, y - 1

        self.axial_scale.configure(from_=self.axial_min, to=self.axial_max)
        self.sagital_scale.configure(from_=self.sagital_min, to=self.sagital_max)
        self.coronal_scale.configure(from_=self.coronal_min, to=self.coronal_max)

        self.axial_var.set((self.axial_min + self.axial_max) // 2)
        self.sagital_var.set((self.sagital_min + self.sagital_max) // 2)
        self.coronal_var.set((self.coronal_min + self.coronal_max) // 2)

        for k in self.zoom_state:
            self.zoom_state[k]["base"] = None
            self.zoom_state[k]["current"] = None

        self.on_axial_change(self.axial_var.get())
        self.on_sagital_change(self.sagital_var.get())
        self.on_coronal_change(self.coronal_var.get())

    # =========================
    # RENDER
    # =========================
    def _get_current_data(self):
        if self.current_key is None or self.current_key not in self.volumes:
            return None
        return self.volumes[self.current_key]["data"]

    def _normalize_for_display(self, img2d: np.ndarray) -> np.ndarray:
        img = np.asarray(img2d, dtype=np.float32)
        lo = np.percentile(img, 1.0)
        hi = np.percentile(img, 99.0)
        if hi <= lo:
            return img
        img = np.clip(img, lo, hi)
        return (img - lo) / (hi - lo)

    def _apply_zoom_after_draw(self, view_key: str, ax):
        state = self.zoom_state[view_key]
        if state["base"] is None:
            base = (ax.get_xlim(), ax.get_ylim())
            state["base"] = base
            state["current"] = base
        else:
            cur = state["current"]
            if cur is not None:
                ax.set_xlim(cur[0])
                ax.set_ylim(cur[1])

    def on_axial_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return
        idx = max(self.axial_min, min(self.axial_max, int(idx)))
        sl = self._normalize_for_display(data[:, :, idx])

        self.axial_axes.clear()
        self.axial_axes.set_axis_off()
        self.axial_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self._apply_zoom_after_draw("axial", self.axial_axes)
        self.axial_canvas.draw_idle()

    def on_sagital_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return
        idx = max(self.sagital_min, min(self.sagital_max, int(idx)))
        sl = self._normalize_for_display(data[idx, :, :])

        self.sagital_axes.clear()
        self.sagital_axes.set_axis_off()
        self.sagital_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self._apply_zoom_after_draw("sagital", self.sagital_axes)
        self.sagital_canvas.draw_idle()

    def on_coronal_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return
        idx = max(self.coronal_min, min(self.coronal_max, int(idx)))
        sl = self._normalize_for_display(data[:, idx, :])

        self.coronal_axes.clear()
        self.coronal_axes.set_axis_off()
        self.coronal_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self._apply_zoom_after_draw("coronal", self.coronal_axes)
        self.coronal_canvas.draw_idle()

    def _render_placeholders(self):
        for key, (ax, canvas) in {
            "axial": (self.axial_axes, self.axial_canvas),
            "sagital": (self.sagital_axes, self.sagital_canvas),
            "coronal": (self.coronal_axes, self.coronal_canvas),
        }.items():
            ax.clear()
            ax.set_axis_off()
            ax.text(0.5, 0.5, "Import NIfTI dulu", ha="center", va="center")
            canvas.draw()
            self.zoom_state[key]["base"] = None
            self.zoom_state[key]["current"] = None

    # =========================
    # SLIDER HELPERS
    # =========================
    def _set_from_slider(self, var: tk.IntVar, val, on_change):
        v = int(float(val))
        if var.get() != v:
            var.set(v)
        on_change(v)

    def _step(self, var: tk.IntVar, delta: int, on_change):
        v = var.get() + delta
        if var is self.axial_var:
            v = max(self.axial_min, min(self.axial_max, v))
        elif var is self.sagital_var:
            v = max(self.sagital_min, min(self.sagital_max, v))
        elif var is self.coronal_var:
            v = max(self.coronal_min, min(self.coronal_max, v))
        var.set(v)
        on_change(v)

    # =========================
    # ZOOM CONTROLS
    # =========================
    def zoom_in(self, view_key: str):
        self._apply_zoom(view_key, factor=0.8)

    def zoom_out(self, view_key: str):
        self._apply_zoom(view_key, factor=1.25)

    def zoom_fit(self, view_key: str):
        ax, canvas = self.view_map[view_key]
        base = self.zoom_state[view_key]["base"]
        if base is None:
            return
        ax.set_xlim(base[0])
        ax.set_ylim(base[1])
        self.zoom_state[view_key]["current"] = base
        canvas.draw_idle()

    def _apply_zoom(self, view_key: str, factor: float):
        ax, canvas = self.view_map[view_key]

        base = self.zoom_state[view_key]["base"]
        if base is None:
            base = (ax.get_xlim(), ax.get_ylim())
            self.zoom_state[view_key]["base"] = base

        current = self.zoom_state[view_key]["current"]
        if current is None:
            current = (ax.get_xlim(), ax.get_ylim())

        (x0, x1), (y0, y1) = current
        cx = (x0 + x1) / 2.0
        cy = (y0 + y1) / 2.0

        w = (x1 - x0) * factor
        h = (y1 - y0) * factor

        new_xlim = (cx - w / 2.0, cx + w / 2.0)
        new_ylim = (cy - h / 2.0, cy + h / 2.0)

        ax.set_xlim(new_xlim)
        ax.set_ylim(new_ylim)

        self.zoom_state[view_key]["current"] = (new_xlim, new_ylim)
        canvas.draw_idle()

    # =========================
    # SAVE PLACEHOLDER
    # =========================
    def on_save(self):
        messagebox.showinfo("Save", "Nanti fungsi save bisa kamu isi untuk export slice atau screenshot.")


if __name__ == "__main__":
    app = CTViewerGUI()
    app.mainloop()
