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

        # Menampung semua volume yang diimport selama aplikasi berjalan
        self.volumes = {}
        self.current_key = None

        # Slider ranges default (akan berubah saat volume dipilih)
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

        self._render_placeholders()

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
    # MAIN LAYOUT (DIUBAH SESUAI SKETSA)
    # =========================
    def _build_main_layout(self):
        # container utama
        container = ttk.Frame(self, style="App.TFrame", padding=(16, 8, 16, 16))
        container.grid(row=1, column=0, sticky="nsew")
        self.grid_rowconfigure(1, weight=1)

        # 2 kolom: kiri tools, kanan visualisasi
        container.grid_columnconfigure(0, weight=0)   # tools
        container.grid_columnconfigure(1, weight=1)   # viewer
        container.grid_rowconfigure(0, weight=1)

        # =================
        # KIRI: TOOL PANEL
        # =================
        tool_panel = ttk.Labelframe(
            container,
            text="AREA TOOL SEGMENTASI",
            style="Panel.TLabelframe",
            padding=(10, 10)
        )
        tool_panel.grid(row=0, column=0, sticky="ns", padx=(0, 12))
        tool_panel.configure(width=230)
        tool_panel.grid_propagate(False)

        # Placeholder isi tool
        ttk.Label(
            tool_panel,
            text="(Placeholder)\n\n• Brush\n• Eraser\n• Threshold\n• Undo/Redo\n• Overlay Mask",
            justify="left"
        ).pack(anchor="nw")

        # ===========================
        # KANAN: VISUALISASI 2x2 GRID
        # ===========================
        vis = ttk.Frame(container, style="App.TFrame")
        vis.grid(row=0, column=1, sticky="nsew")

        vis.grid_columnconfigure(0, weight=1)
        vis.grid_columnconfigure(1, weight=1)
        # baris 0 slider atas (axial, sagital)
        vis.grid_rowconfigure(0, weight=0)
        # baris 1 view atas (axial, sagital)
        vis.grid_rowconfigure(1, weight=1)
        # baris 2 slider bawah (coronal di kiri, kosong di kanan)
        vis.grid_rowconfigure(2, weight=0)
        # baris 3 view bawah (coronal, 3d)
        vis.grid_rowconfigure(3, weight=1)

        # Baris 0: slider axial dan sagital
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

        # Baris 1: view axial dan sagital
        self.axial_panel, self.axial_axes, self.axial_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI AXIAL"
        )
        self.axial_panel.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        self.sagital_panel, self.sagital_axes, self.sagital_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI SAGITAL"
        )
        self.sagital_panel.grid(row=1, column=1, sticky="nsew", padx=4, pady=4)

        # Baris 2: slider coronal (kiri), kanan kosong
        self.coronal_slider_row = self._build_slider_row(
            parent=vis,
            label_text="SLIDER UNTUK GANTI SLICE CORONAL",
            var=self.coronal_var,
            store_scale_as="coronal_scale",
            on_change=self.on_coronal_change
        )
        self.coronal_slider_row.grid(row=2, column=0, sticky="ew", padx=4, pady=(6, 6))

        ttk.Frame(vis, style="App.TFrame").grid(row=2, column=1, sticky="ew")

        # Baris 3: view coronal dan 3D
        self.coronal_panel, self.coronal_axes, self.coronal_canvas = self._build_view_plot_only(
            parent=vis,
            title="GRID UNTUK VISUALISASI CORONAL"
        )
        self.coronal_panel.grid(row=3, column=0, sticky="nsew", padx=4, pady=4)

        panel3d = ttk.Labelframe(vis, text="GRID UNTUK 3D", style="Panel.TLabelframe", padding=(10, 10))
        panel3d.grid(row=3, column=1, sticky="nsew", padx=4, pady=4)
        panel3d.grid_rowconfigure(0, weight=1)
        panel3d.grid_columnconfigure(0, weight=1)
        ttk.Label(panel3d, text="TAMPILAN 3D (placeholder)", anchor="center").grid(row=0, column=0, sticky="nsew")

    # =========================
    # SLIDER ROW (LABEL + ◀ + SLIDER + ▶ + VALUE)
    # =========================
    def _build_slider_row(self, parent, label_text, var, store_scale_as, on_change):
        frame = ttk.Frame(parent, style="App.TFrame")
        frame.grid_columnconfigure(2, weight=1)

        ttk.Label(frame, text=label_text).grid(row=0, column=0, sticky="w", padx=(0, 10))

        btn_left = ttk.Button(
            frame,
            text="◀",
            width=3,
            command=lambda: self._step(var, -1, on_change)
        )
        btn_left.grid(row=0, column=1, sticky="w")

        slider = ttk.Scale(
            frame,
            from_=0,
            to=10,
            orient="horizontal",
            command=lambda val: self._set_from_slider(var, val, on_change)
        )
        slider.grid(row=0, column=2, sticky="ew", padx=(8, 8))

        btn_right = ttk.Button(
            frame,
            text="▶",
            width=3,
            command=lambda: self._step(var, 1, on_change)
        )
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
    # VIEW ONLY (PANEL + MATPLOTLIB CANVAS)
    # =========================
    def _build_view_plot_only(self, parent, title):
        panel = ttk.Labelframe(parent, text=title, style="Panel.TLabelframe", padding=(6, 6))
        panel.grid_rowconfigure(0, weight=1)
        panel.grid_columnconfigure(0, weight=1)

        fig = Figure(figsize=(4, 4), dpi=100)
        ax = fig.add_subplot(111)
        ax.set_axis_off()

        canvas = FigureCanvasTkAgg(fig, master=panel)
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        return panel, ax, canvas

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

        base = os.path.basename(path)
        key = self._make_unique_key(base)

        self.volumes[key] = {
            "path": path,
            "data": data,
            "shape": data.shape
        }

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

        x, y, z = data.shape

        # Axial sepanjang Z, Sagital sepanjang X, Coronal sepanjang Y
        self.axial_min, self.axial_max = 0, z - 1
        self.sagital_min, self.sagital_max = 0, x - 1
        self.coronal_min, self.coronal_max = 0, y - 1

        self.axial_scale.configure(from_=self.axial_min, to=self.axial_max)
        self.sagital_scale.configure(from_=self.sagital_min, to=self.sagital_max)
        self.coronal_scale.configure(from_=self.coronal_min, to=self.coronal_max)

        self.axial_var.set((self.axial_min + self.axial_max) // 2)
        self.sagital_var.set((self.sagital_min + self.sagital_max) // 2)
        self.coronal_var.set((self.coronal_min + self.coronal_max) // 2)

        self.on_axial_change(self.axial_var.get())
        self.on_sagital_change(self.sagital_var.get())
        self.on_coronal_change(self.coronal_var.get())

    # =========================
    # RENDER
    # =========================
    def _get_current_data(self):
        if self.current_key is None:
            return None
        if self.current_key not in self.volumes:
            return None
        return self.volumes[self.current_key]["data"]

    def _normalize_for_display(self, img2d: np.ndarray) -> np.ndarray:
        img = np.asarray(img2d, dtype=np.float32)
        lo = np.percentile(img, 1.0)
        hi = np.percentile(img, 99.0)
        if hi <= lo:
            return img
        img = np.clip(img, lo, hi)
        img = (img - lo) / (hi - lo)
        return img

    def on_axial_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return

        idx = int(idx)
        idx = max(self.axial_min, min(self.axial_max, idx))

        sl = data[:, :, idx]
        sl = self._normalize_for_display(sl)

        self.axial_axes.clear()
        self.axial_axes.set_axis_off()
        self.axial_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self.axial_canvas.draw_idle()

    def on_sagital_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return

        idx = int(idx)
        idx = max(self.sagital_min, min(self.sagital_max, idx))

        sl = data[idx, :, :]
        sl = self._normalize_for_display(sl)

        self.sagital_axes.clear()
        self.sagital_axes.set_axis_off()
        self.sagital_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self.sagital_canvas.draw_idle()

    def on_coronal_change(self, idx: int):
        data = self._get_current_data()
        if data is None:
            return

        idx = int(idx)
        idx = max(self.coronal_min, min(self.coronal_max, idx))

        sl = data[:, idx, :]
        sl = self._normalize_for_display(sl)

        self.coronal_axes.clear()
        self.coronal_axes.set_axis_off()
        self.coronal_axes.imshow(np.rot90(sl), cmap="gray", interpolation="nearest")
        self.coronal_canvas.draw_idle()

    def _render_placeholders(self):
        self.axial_axes.clear()
        self.axial_axes.set_axis_off()
        self.axial_axes.text(0.5, 0.5, "Import NIfTI dulu", ha="center", va="center")
        self.axial_canvas.draw()

        self.sagital_axes.clear()
        self.sagital_axes.set_axis_off()
        self.sagital_axes.text(0.5, 0.5, "Import NIfTI dulu", ha="center", va="center")
        self.sagital_canvas.draw()

        self.coronal_axes.clear()
        self.coronal_axes.set_axis_off()
        self.coronal_axes.text(0.5, 0.5, "Import NIfTI dulu", ha="center", va="center")
        self.coronal_canvas.draw()

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
    # SAVE PLACEHOLDER
    # =========================
    def on_save(self):
        messagebox.showinfo("Save", "Nanti fungsi save bisa kamu isi untuk export slice atau screenshot.")


if __name__ == "__main__":
    app = CTViewerGUI()
    app.mainloop()
