from sod1_analysis import analyze_fasta
import sys
import json
import subprocess
from pathlib import Path

from PySide6.QtCore import Qt, QSettings
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFrame, QLabel,
    QPushButton, QVBoxLayout, QHBoxLayout, QGridLayout,
    QStackedWidget, QScrollArea, QTableWidget,
    QTableWidgetItem, QHeaderView, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QMessageBox
)

PROJECT = "/home/kyokugen/protein-work/ProteinMPNN/sod1_project"
MANIFEST = f"{PROJECT}/results/test10/manifest.json"
FASTA = f"{PROJECT}/results/test10/seqs/2C9V.fa"


def read_wsl_file(path):
    result = subprocess.run(
        ["wsl.exe", "cat", path],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=20,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())
    return result.stdout


def load_research():
    manifest = json.loads(read_wsl_file(MANIFEST))
    fasta = read_wsl_file(FASTA)

    records = []
    sequence = ""
    for line in fasta.splitlines():
        if line.startswith(">"):
            if sequence:
                records.append(sequence)
            sequence = ""
        else:
            sequence += line.strip()
    if sequence:
        records.append(sequence)

    samples = records[1:]
    fixed = {
        46: "H", 48: "H", 57: "C", 63: "H",
        71: "H", 80: "H", 83: "D", 120: "H",
        124: "D", 143: "R", 146: "C"
    }

    def valid(seq):
        chains = seq.split("/")
        return (
            len(chains) == 2
            and all(len(s) == 153 for s in chains)
            and chains[0] == chains[1]
            and all(
                s[pos - 1] == aa
                for s in chains
                for pos, aa in fixed.items()
            )
        )

    return {
        "manifest": manifest,
        "samples": samples,
        "valid": sum(valid(s) for s in samples),
        "unique": len(set(s.split("/")[0] for s in samples)),
        "validator": valid
    }


def label(text, size=13, bold=False, color="#18233B"):
    widget = QLabel(str(text))
    widget.setWordWrap(True)
    widget.setStyleSheet(
        f"font-size:{size}px;"
        f"font-weight:{700 if bold else 400};"
        f"color:{color};"
    )
    return widget


def card():
    frame = QFrame()
    frame.setObjectName("card")
    return frame


TRANSLATIONS = {'研究ダッシュボード': 'Research Dashboard', 'タンパク質の配列設計と解析を、ひとつの場所で。': 'Protein sequence design and analysis in one place.', 'プロジェクト': 'Project', 'SOD1研究': 'SOD1 Research', '構造（PDB）': 'Structure (PDB)', '生成した配列': 'Designed Sequences', '検証済み': 'Validated', 'まずはここから': 'Getting Started', '研究の流れに沿って進められます。': 'Follow the research workflow step by step.', '① 構造を読み込む': '① Load Structure', '② 配列を設計する': '② Design Sequences', '③ 結果を確認する': '③ Review Results', '④ 研究記録を残す': '④ Save Records', 'PDBファイル': 'PDB File', '配列・制約検証': 'Sequences / Validation', '条件・結果': 'Conditions / Results', '保存済みの試験結果': 'Saved Test Results', '配列設計に使用するタンパク質構造を確認します。': 'Review the protein structure used for sequence design.', '現在の構造': 'Current Structure', '2C9V / ヒトSOD1': '2C9V / Human SOD1', '対象鎖：A鎖・F鎖（各153残基）': 'Target chains: A and F (153 residues each)', '※ 現在は保存済み構造の情報を表示しています。': 'Note: Showing information from the saved structure.', '配列設計': 'Sequence Design', 'ProteinMPNNの生成条件を設定します。': 'Configure ProteinMPNN generation parameters.', '基本設定': 'Basic Settings', '生成する配列数': 'Number of Sequences', 'Temperature（多様性）': 'Temperature (Diversity)', 'Seed（乱数の初期値）': 'Seed (Random Initialization)', '設計制約': 'Design Constraints', '対象：ヒトSOD1 / 2C9V': 'Target: Human SOD1 / 2C9V', '設計鎖：A鎖・F鎖': 'Design chains: A and F', '二量体制約：153組': 'Dimer constraints: 153 pairs', '固定残基：11か所': 'Fixed residues: 11 positions', '待機中': 'Ready', '配列設計を開始': 'Start Sequence Design', '設定内容を確認する': 'Review Settings', '研究記録': 'Research Records', '再現性を確保するための実行環境と条件を確認します。': 'Review the environment and parameters for reproducibility.', '研究記録を読み込めませんでした。': 'Could not load research records.', '不明': 'Unknown', '計算に失敗しました': 'Computation failed', 'ProteinMPNNを実行中です。完了までお待ちください。': 'ProteinMPNN is running. Please wait.'}

class Studio(QMainWindow):
    def __init__(self):
        super().__init__()
        self.settings = QSettings(
            "ProteinMPNN", "ResearchStudio"
        )
        self.language = self.settings.value(
            "language", "ja", type=str
        )
        if self.language not in ("ja", "en"):
            self.language = "ja"
        self.setWindowTitle("ProteinMPNN Research Studio")
        self.resize(1250, 800)

        try:
            self.data = load_research()
            self.error = None
        except Exception as exc:
            self.data = None
            self.error = str(exc)

        root = QWidget()
        self.setCentralWidget(root)
        main = QHBoxLayout(root)
        main.setContentsMargins(0, 0, 0, 0)
        main.setSpacing(0)

        sidebar = QFrame()
        sidebar.setObjectName("sidebar")
        sidebar.setFixedWidth(255)
        side = QVBoxLayout(sidebar)
        side.setContentsMargins(20, 30, 20, 22)
        side.setSpacing(12)

        side.addWidget(label("🧬  ProteinMPNN", 19, True))
        side.addWidget(label("Research Studio", 14, True))
        side.addSpacing(28)

        self.pages = QStackedWidget()
        self.buttons = []

        self.menu_names = [
            ("⌂   ホーム", "⌂   Home"),
            ("◇   構造（PDB）", "◇   Structure (PDB)"),
            ("⚙   配列設計", "⚙   Sequence Design"),
            ("▤   結果の確認・解析", "▤   Results"),
            ("▣   研究記録", "▣   Research Records"),
            ("⚙   設定", "⚙   Settings"),
        ]
        builders = [
            self.home_page,
            self.structure_page,
            self.design_page,
            self.results_page,
            self.records_page,
            self.settings_page,
        ]
        menus = [
            (names[0 if self.language == "ja" else 1], builder)
            for names, builder in zip(self.menu_names, builders)
        ]

        for index, (name, builder) in enumerate(menus):
            button = QPushButton(name)
            button.setObjectName("nav")
            button.setCheckable(True)
            button.setMinimumHeight(46)
            button.clicked.connect(
                lambda checked=False, i=index: self.navigate(i)
            )
            self.buttons.append(button)
            side.addWidget(button)
            self.pages.addWidget(self.scroll_page(builder()))

        side.addStretch()

        self.language_label = QLabel(
            "表示言語" if self.language == "ja"
            else "Language"
        )
        self.language_label.setStyleSheet(
            "color: #69788F; font-size: 12px;"
        )
        side.addWidget(self.language_label)

        self.language_selector = QComboBox()
        self.language_selector.addItem("日本語", "ja")
        self.language_selector.addItem("English", "en")
        self.language_selector.setCurrentIndex(
            0 if self.language == "ja" else 1
        )
        self.language_selector.setMinimumHeight(38)
        side.addWidget(self.language_selector)

        status = card()
        status_layout = QVBoxLayout(status)
        self.engine_label = label(
            "計算エンジン" if self.language == "ja"
            else "Compute Engine", 13, True
        )
        status_layout.addWidget(self.engine_label)
        status_layout.addWidget(label("WSL2 / Ubuntu", 12))
        self.engine_status = label(
            "", 12, False,
            "#16834A" if self.data else "#B45309"
        )
        status_layout.addWidget(self.engine_status)
        self.update_language_labels()
        side.addWidget(status)

        self.language_selector.currentIndexChanged.connect(
            self.change_language
        )

        main.addWidget(sidebar)
        main.addWidget(self.pages, 1)

        self.setStyleSheet("""
            QWidget {
                background: #F7F9FC;
                color: #18233B;
                font-family: "Yu Gothic UI", "Meiryo";
                font-size: 13px;
            }
            QFrame#sidebar {
                background: white;
                border-right: 1px solid #E5EAF2;
            }
            QFrame#card {
                background: white;
                border: 1px solid #E5EAF2;
                border-radius: 12px;
            }
            QFrame#card QLabel {
                background: transparent;
                border: none;
            }
            QPushButton#nav {
                text-align: left;
                padding: 10px 15px;
                background: transparent;
                border: none;
                border-radius: 9px;
                font-size: 13px;
            }
            QPushButton#nav:checked {
                background: #E5EEFF;
                color: #175CE5;
                font-weight: bold;
            }
            QPushButton#nav:hover {
                background: #EDF3FF;
            }
            QPushButton#action {
                background: #2468E8;
                color: white;
                border: none;
                border-radius: 9px;
                padding: 12px 18px;
                font-weight: bold;
            }
            QPushButton#action:disabled {
                background: #E7ECF3;
                color: #8793A5;
            }
            QTableWidget {
                background: white;
                border: 1px solid #E5EAF2;
                border-radius: 9px;
                gridline-color: #EDF1F6;
            }
        """)

        self.navigate(0)

    def _latest_design_fasta(self):
        root = PROJECT + "/results"
        command = [
            "wsl.exe", "-d", "Ubuntu", "--", "find", root,
            "-maxdepth", "3", "-type", "f", "-name", "2C9V.fa"
        ]
        process = subprocess.run(
            command, capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=30,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)
        )
        if process.returncode != 0:
            raise RuntimeError(process.stderr.strip() or "FASTA search failed")
        paths = [p.strip() for p in process.stdout.splitlines() if p.strip()]
        gui_paths = [p for p in paths if "/results/gui_" in p]
        if gui_paths:
            return max(gui_paths)
        if FASTA in paths:
            return FASTA
        raise FileNotFoundError("2C9V.fa was not found in results")

    def _refresh_results_page(self):
        index = 3
        old = self.pages.widget(index)
        new = self.scroll_page(self.results_page())
        self.pages.insertWidget(index, new)
        if old is not None:
            self.pages.removeWidget(old)
            old.deleteLater()
        self.pages.setCurrentIndex(index)
        for i, button in enumerate(self.buttons):
            button.setChecked(i == index)

    def tr_ui(self, ja, en):
        return ja if self.language == "ja" else en

    def update_language_labels(self):
        for button, names in zip(self.buttons, self.menu_names):
            button.setText(
                names[0 if self.language == "ja" else 1]
            )

        self.language_label.setText(
            self.tr_ui("表示言語", "Language")
        )
        self.engine_label.setText(
            self.tr_ui("計算エンジン", "Compute Engine")
        )
        if self.data:
            self.engine_status.setText(
                self.tr_ui(
                    "● 研究データ読み込みOK",
                    "● Research data loaded"
                )
            )
        else:
            self.engine_status.setText(
                self.tr_ui(
                    "● 接続・読み込みを確認",
                    "● Check connection / data"
                )
            )

    def change_language(self, index):
        language = self.language_selector.itemData(index)
        if language not in ("ja", "en"):
            return
        if language == self.language:
            return

        self.language = language
        self.settings.setValue("language", language)
        self.update_language_labels()
        for page_index in (0, 1, 2, 4):
            self.translate_page_widgets(self.pages.widget(page_index))

        current = self.pages.currentIndex()
        if current == 3:
            self._refresh_results_page()
        elif current == 5:
            old = self.pages.widget(5)
            new = self.scroll_page(self.settings_page())
            self.pages.insertWidget(5, new)
            self.pages.removeWidget(old)
            old.deleteLater()
            self.pages.setCurrentIndex(5)

    def translate_display_text(self, text):
        if self.language != "en":
            return text
        translations = TRANSLATIONS
        if text in translations:
            return translations[text]
        if text.startswith("固定残基："):
            return "Fixed residues: " + text[len("固定残基："):]
        if text.startswith("記録ファイル数: "):
            return "Recorded files: " + text[len("記録ファイル数: "):]
        if text.startswith("研究データを読み込めませんでした。\n"):
            return "Could not load research data.\n" + text.split("\n", 1)[1]
        if text.startswith("2C9V  |  生成 "):
            return (text.replace("生成 ", "Generated ")
                        .replace("検証合格 ", "Validated ")
                        .replace("重複を除く ", "Unique ")
                        .replace(" 配列", " sequences"))
        return text

    def translate_page_widgets(self, page):
        for widget in [page] + page.findChildren(QWidget):
            if isinstance(widget, QLabel):
                original = widget.property("i18n_original")
                if original is None:
                    original = widget.text()
                    widget.setProperty("i18n_original", original)
                widget.setText(self.translate_display_text(original))
            elif isinstance(widget, QPushButton):
                # Sidebar buttons are translated separately.
                if widget in self.buttons:
                    continue
                original = widget.property("i18n_original")
                if original is None:
                    original = widget.text()
                    widget.setProperty("i18n_original", original)
                widget.setText(self.translate_display_text(original))

    def navigate(self, index):
        if index == 3:
            self._refresh_results_page()
            return
        self.pages.setCurrentIndex(index)
        for i, button in enumerate(self.buttons):
            button.setChecked(i == index)

    def scroll_page(self, content):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidget(content)
        self.translate_page_widgets(content)
        return scroll

    def page(self, title, subtitle):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(32, 32, 32, 32)
        layout.setSpacing(22)
        layout.addWidget(label(title, 27, True))
        layout.addWidget(label(subtitle, 13, False, "#69788F"))
        return widget, layout

    def create_metric_card(self, title, value):
        frame = card()
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.addWidget(label(title, 12, False, "#69788F"))
        layout.addWidget(label(value, 25, True))
        return frame

    def home_page(self):
        widget, layout = self.page(
            "研究ダッシュボード",
            "タンパク質の配列設計と解析を、ひとつの場所で。"
        )

        row = QHBoxLayout()
        values = [
            ("プロジェクト", "SOD1研究"),
            ("構造（PDB）", "2C9V"),
            ("生成した配列", str(len(self.data["samples"])) if self.data else "—"),
            ("検証済み", str(self.data["valid"]) if self.data else "—")
        ]
        for title, value in values:
            row.addWidget(self.create_metric_card(title, value))
        layout.addLayout(row)

        section = card()
        inner = QVBoxLayout(section)
        inner.setContentsMargins(24, 22, 24, 24)
        inner.setSpacing(18)
        inner.addWidget(label("まずはここから", 19, True))
        inner.addWidget(label("研究の流れに沿って進められます。", 12, False, "#69788F"))

        steps = QHBoxLayout()
        for title, description in [
            ("① 構造を読み込む", "PDBファイル"),
            ("② 配列を設計する", "ProteinMPNN"),
            ("③ 結果を確認する", "配列・制約検証"),
            ("④ 研究記録を残す", "条件・結果")
        ]:
            step = card()
            sl = QVBoxLayout(step)
            sl.setContentsMargins(12, 18, 12, 18)
            sl.addWidget(label(title, 13, True))
            sl.addWidget(label(description, 11, False, "#69788F"))
            steps.addWidget(step)
        inner.addLayout(steps)
        layout.addWidget(section)

        recent = card()
        rl = QVBoxLayout(recent)
        rl.setContentsMargins(24, 22, 24, 22)
        rl.addWidget(label("保存済みの試験結果", 18, True))
        if self.data:
            rl.addWidget(label(
                f"2C9V  |  生成 {len(self.data['samples'])} 配列"
                f"  |  検証合格 {self.data['valid']} 配列"
                f"  |  重複を除く {self.data['unique']} 配列"
            ))
        else:
            rl.addWidget(label(
                f"研究データを読み込めませんでした。\n{self.error}",
                12, False, "#B45309"
            ))
        layout.addWidget(recent)
        layout.addStretch()
        return widget

    def structure_page(self):
        widget, layout = self.page(
            "構造（PDB）",
            "配列設計に使用するタンパク質構造を確認します。"
        )
        layout.addWidget(self.create_metric_card("現在の構造", "2C9V / ヒトSOD1"))
        layout.addWidget(label("対象鎖：A鎖・F鎖（各153残基）"))
        layout.addWidget(label(
            "固定残基：H46, H48, C57, H63, H71, H80, "
            "D83, H120, D124, R143, C146"
        ))
        layout.addWidget(label(
            "※ 現在は保存済み構造の情報を表示しています。",
            12, False, "#69788F"
        ))
        layout.addStretch()
        return widget


    def design_page(self):
        widget, layout = self.page(
            "配列設計",
            "ProteinMPNNの生成条件を設定します。"
        )

        parameters = (
            self.data["manifest"].get("test_parameters", {})
            if self.data else {}
        )

        form_card = card()
        form = QVBoxLayout(form_card)
        form.setContentsMargins(28, 24, 28, 24)
        form.setSpacing(14)

        form.addWidget(label("基本設定", 19, True))

        form.addWidget(label("生成する配列数", 14, True))
        self.num_sequences = QSpinBox()
        self.num_sequences.setRange(1, 100000)
        self.num_sequences.setValue(
            int(parameters.get("num_seq_per_target", 10))
        )
        self.num_sequences.setMinimumHeight(42)
        form.addWidget(self.num_sequences)

        form.addWidget(label("Temperature（多様性）", 14, True))
        self.temperature = QDoubleSpinBox()
        self.temperature.setRange(0.01, 1.0)
        self.temperature.setDecimals(2)
        self.temperature.setSingleStep(0.05)
        self.temperature.setValue(
            float(parameters.get("sampling_temp", 0.1))
        )
        self.temperature.setMinimumHeight(42)
        form.addWidget(self.temperature)

        form.addWidget(label("Seed（乱数の初期値）", 14, True))
        self.seed_input = QSpinBox()
        self.seed_input.setRange(1, 2147483647)
        self.seed_input.setValue(
            int(parameters.get("seed", 42))
        )
        self.seed_input.setMinimumHeight(42)
        form.addWidget(self.seed_input)

        # Model configuration (SOD1 preset; retains existing constraints).
        form.addWidget(label(self.tr_ui("モデルの種類", "Model Type"), 14, True))
        self.model_type_combo = QComboBox()
        self.model_type_combo.addItem("Vanilla ProteinMPNN", "vanilla")
        self.model_type_combo.addItem("Soluble ProteinMPNN", "soluble")
        self.model_type_combo.addItem("Cα-only ProteinMPNN", "ca")
        self.model_type_combo.setMinimumHeight(42)
        form.addWidget(self.model_type_combo)

        form.addWidget(label(self.tr_ui("モデルの重み", "Model Weights"), 14, True))
        self.model_weight_combo = QComboBox()
        self.model_weight_combo.setMinimumHeight(42)
        form.addWidget(self.model_weight_combo)

        def refresh_model_weights():
            kind = self.model_type_combo.currentData()
            previous = self.model_weight_combo.currentText()
            weights = ["v_48_002", "v_48_010", "v_48_020"]
            if kind != "ca":
                weights.append("v_48_030")
            self.model_weight_combo.clear()
            self.model_weight_combo.addItems(weights)
            chosen = previous if previous in weights else "v_48_020"
            self.model_weight_combo.setCurrentText(chosen)

        self.model_type_combo.currentIndexChanged.connect(refresh_model_weights)
        refresh_model_weights()

        form.addWidget(label("Batch Size", 14, True))
        self.batch_size_input = QSpinBox()
        self.batch_size_input.setRange(1, 100000)
        self.batch_size_input.setValue(1)
        self.batch_size_input.setMinimumHeight(42)
        form.addWidget(self.batch_size_input)

        form.addWidget(label("Backbone Noise", 14, True))
        self.backbone_noise_input = QDoubleSpinBox()
        self.backbone_noise_input.setRange(0.0, 2.0)
        self.backbone_noise_input.setDecimals(2)
        self.backbone_noise_input.setSingleStep(0.05)
        self.backbone_noise_input.setValue(0.0)
        self.backbone_noise_input.setMinimumHeight(42)
        form.addWidget(self.backbone_noise_input)

        layout.addWidget(form_card)

        constraint_card = card()
        constraints = QVBoxLayout(constraint_card)
        constraints.setContentsMargins(28, 24, 28, 24)

        constraints.addWidget(label("設計制約", 19, True))
        constraints.addWidget(label("対象：ヒトSOD1 / 2C9V"))
        constraints.addWidget(label("設計鎖：A鎖・F鎖"))
        constraints.addWidget(label("二量体制約：153組"))
        constraints.addWidget(label("固定残基：11か所"))
        constraints.addWidget(label(
            "H46 / H48 / C57 / H63 / H71 / H80 / "
            "D83 / H120 / D124 / R143 / C146"
        ))

        layout.addWidget(constraint_card)

        self.run_status = QLabel("待機中")
        layout.addWidget(self.run_status)
        self.run_button = QPushButton("配列設計を開始")
        self.run_button.setObjectName("action")
        self.run_button.setMinimumHeight(48)
        self.run_button.clicked.connect(self.start_design_run)
        layout.addWidget(self.run_button)

        button = QPushButton("設定内容を確認する")
        button.setObjectName("action")
        button.setMinimumHeight(48)
        button.clicked.connect(self.confirm_design)
        layout.addWidget(button)

        layout.addStretch()
        return widget

    def confirm_design(self):
        message = "\n".join([
            f"生成配列数：{self.num_sequences.value()}",
            f"Temperature：{self.temperature.value():.2f}",
            f"Seed：{self.seed_input.value()}",
            f"Model：{self.model_type_combo.currentText()} / {self.model_weight_combo.currentText()}",
            f"Batch Size：{self.batch_size_input.value()}",
            f"Backbone Noise：{self.backbone_noise_input.value():.2f}",
            "",
            "現在は設定確認のみです。",
            "配列生成の計算はまだ実行されません。"
        ])

        QMessageBox.information(
            self,
            "配列設計の設定確認",
            message
        )

    def start_design_run(self):
        from mpnn_runner_extended import run_mpnn
        from PySide6.QtCore import QThread, Signal

        count = self.num_sequences.value()
        temp = self.temperature.value()
        seed = self.seed_input.value()
        model_type = self.model_type_combo.currentData()
        model_name = self.model_weight_combo.currentText()
        batch_size = self.batch_size_input.value()
        backbone_noise = self.backbone_noise_input.value()
        if count % batch_size != 0:
            QMessageBox.warning(
                self,
                self.tr_ui("設定エラー", "Invalid settings"),
                self.tr_ui(
                    "生成配列数はBatch Sizeで割り切れる値にしてください。",
                    "The number of sequences must be divisible by Batch Size."
                )
            )
            return
        answer = QMessageBox.question(
            self, "配列設計の実行確認",
            f"{count}配列を生成します。\nTemperature: {temp:.2f}\nSeed: {seed}\n実行しますか？"
        )
        if answer != QMessageBox.StandardButton.Yes:
            return

        class Worker(QThread):
            finished_ok = Signal(str)
            failed = Signal(str)

            def run(self):
                try:
                    output = run_mpnn(
                        count, temp, seed,
                        model_name=model_name,
                        batch_size=batch_size,
                        backbone_noise=backbone_noise,
                        use_soluble_model=(model_type == "soluble"),
                        ca_only=(model_type == "ca"),
                    )
                    fasta = read_wsl_file(output + "/seqs/2C9V.fa")
                    fixed = {46:"H",48:"H",57:"C",63:"H",71:"H",80:"H",83:"D",120:"H",124:"D",143:"R",146:"C"}
                    records = []
                    header = None
                    chunks = []
                    for line in fasta.splitlines():
                        if line.startswith(">"):
                            if header is not None:
                                records.append((header, "".join(chunks)))
                            header, chunks = line, []
                        else:
                            chunks.append(line.strip())
                    if header is not None:
                        records.append((header, "".join(chunks)))
                    samples = [seq for h, seq in records if h.startswith(">T=")]
                    valid = 0
                    for seq in samples:
                        chains = seq.split("/")
                        if (len(chains) == 2 and all(len(c) == 153 for c in chains)
                            and chains[0] == chains[1]
                            and all(c[pos-1] == aa for c in chains for pos, aa in fixed.items())):
                            valid += 1
                    if len(samples) != count or valid != count:
                        raise RuntimeError(f"FASTA validation failed: {valid}/{len(samples)} PASS; expected {count}. Output: {output}")
                    self.finished_ok.emit(output + f"\n検証: {valid}/{len(samples)} PASS")
                except Exception as exc:
                    self.failed.emit(str(exc))

        self.run_button.setEnabled(False)
        self.run_status.setText(self.translate_display_text("ProteinMPNNを実行中です。完了までお待ちください。"))
        self.mpnn_worker = Worker(self)
        self.mpnn_worker.finished_ok.connect(self.on_design_success)
        self.mpnn_worker.failed.connect(self.on_design_failure)
        self.mpnn_worker.start()

    def on_design_success(self, output):
        self.run_button.setEnabled(True)
        self.run_status.setText("計算完了: " + output)
        try:
            self._refresh_results_page()
            QMessageBox.information(
                self, "計算完了", "FASTAの生成と解析が完了しました。\n" + output
            )
        except Exception as exc:
            QMessageBox.warning(
                self, "計算完了・解析画面更新失敗",
                "FASTAは生成されましたが、解析画面の更新に失敗しました。\n"
                + str(exc) + "\n" + output
            )

    def on_design_failure(self, error):
        self.run_button.setEnabled(True)
        self.run_status.setText(self.translate_display_text("計算に失敗しました"))
        QMessageBox.critical(self, "実行エラー", error)

    def results_page(self):
        widget, layout = self.page(
            self.tr_ui("配列解析", "Sequence Explorer"),
            self.tr_ui(
                "SOD1設計配列の比較・解析",
                "Explore and compare designed SOD1 variants."
            )
        )

        try:
            fasta_path = self._latest_design_fasta()
            data = analyze_fasta(read_wsl_file(fasta_path))
        except Exception as exc:
            layout.addWidget(label(
                "解析データを読み込めませんでした: " + str(exc)
            ))
            layout.addStretch()
            return widget

        # Summary cards
        stats = QHBoxLayout()
        stats.setSpacing(12)

        for heading, value in [
            (self.tr_ui("生成配列数", "DESIGNED"), data["total"]),
            (self.tr_ui("制約PASS", "VALIDATED"), data["valid"]),
            (self.tr_ui("ユニーク配列", "UNIQUE"), data["unique"]),
        ]:
            card = QFrame()
            card.setObjectName("resultStatCard")
            card.setStyleSheet("""
                QFrame#resultStatCard {
                    background: #F4F7FD;
                    border: 1px solid #E6EBF5;
                    border-radius: 14px;
                }
            """)

            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(18, 15, 18, 15)
            card_layout.setSpacing(5)

            heading_label = QLabel(heading)
            heading_label.setStyleSheet(
                "color: #8390A7; font-size: 11px;"
                "font-weight: 600; border: none;"
            )

            value_label = QLabel(str(value))
            value_label.setStyleSheet(
                "color: #25344D; font-size: 26px;"
                "font-weight: 700; border: none;"
            )

            card_layout.addWidget(heading_label)
            card_layout.addWidget(value_label)
            stats.addWidget(card)

        layout.addLayout(stats)

        # Search and toolbar
        toolbar = QHBoxLayout()
        toolbar.setSpacing(10)

        search_box = QLineEdit()
        search_box.setPlaceholderText(
            self.tr_ui(
                "配列番号・変異位置を検索...",
                "Search variants, positions, mutations..."
            )
        )
        search_box.setClearButtonEnabled(True)
        search_box.setMinimumHeight(42)
        search_box.setStyleSheet("""
            QLineEdit {
                background: white;
                border: 1px solid #DFE6F2;
                border-radius: 11px;
                padding: 8px 14px;
                color: #25344D;
                font-size: 13px;
            }
            QLineEdit:focus {
                border: 1px solid #91A9E8;
            }
        """)

        toolbar.addWidget(search_box, 1)

        count_label = QLabel()
        count_label.setStyleSheet(
            "color: #8390A7; font-size: 12px;"
        )
        toolbar.addWidget(count_label)

        layout.addLayout(toolbar)

        # Results table
        table = QTableWidget(len(data["results"]), 8)
        table.setSortingEnabled(False)
        table.setHorizontalHeaderLabels([
            *(
                [
                    "配列番号", "制約", "変異数", "Score",
                    "Global score", "Recovery", "重複数",
                    "変異位置"
                ]
                if self.language == "ja"
                else [
                    "Variant", "Status", "Mutations", "Score",
                    "Global score", "Recovery", "Duplicates",
                    "Mutation positions"
                ]
            )
        ])

        table.setAlternatingRowColors(True)
        table.setShowGrid(False)
        table.verticalHeader().setVisible(False)
        table.verticalHeader().setDefaultSectionSize(39)
        table.setSelectionBehavior(QTableWidget.SelectRows)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.setMinimumHeight(460)

        table.setStyleSheet("""
            QTableWidget {
                background: white;
                alternate-background-color: #F8FAFE;
                border: 1px solid #E6EBF5;
                border-radius: 12px;
                color: #334155;
                selection-background-color: #E8EEFF;
                selection-color: #25344D;
            }
            QHeaderView::section {
                background: #F1F5FB;
                color: #64748B;
                border: none;
                border-bottom: 1px solid #E6EBF5;
                padding: 10px 8px;
                font-weight: 600;
            }
        """)

        for row, record in enumerate(data["results"]):
            def fmt(value):
                return "-" if value is None else f"{value:.4f}"

            recovery = record["recovery"]
            values = [
                str(record["number"]),
                "PASS" if record["valid"] else "FAIL",
                str(record["mutation_count"]),
                fmt(record["score"]),
                fmt(record["global_score"]),
                "-" if recovery is None
                else f"{recovery * 100:.2f}%",
                str(record["duplicate_count"]),
                record["mutations"],
            ]

            for column, value in enumerate(values):
                item = QTableWidgetItem(value)

                numeric = {
                    0: record["number"],
                    2: record["mutation_count"],
                    3: record["score"],
                    4: record["global_score"],
                    5: None if recovery is None
                       else recovery * 100,
                    6: record["duplicate_count"],
                }

                if column in numeric and numeric[column] is not None:
                    item.setData(Qt.EditRole, float(numeric[column]))

                table.setItem(row, column, item)

        table.horizontalHeader().setSectionResizeMode(
            7, QHeaderView.Stretch
        )
        table.setSortingEnabled(True)

        def filter_results(query):
            query = query.strip().lower()
            visible = 0

            for row in range(table.rowCount()):
                searchable = " ".join(
                    table.item(row, col).text().lower()
                    for col in (0, 1, 7)
                    if table.item(row, col) is not None
                )

                matches = query in searchable
                table.setRowHidden(row, not matches)

                if matches:
                    visible += 1

            count_label.setText(
                self.tr_ui(
                    f"{visible:,} / {table.rowCount():,} 配列",
                    f"{visible:,} / {table.rowCount():,} variants"
                )
            )

        search_box.textChanged.connect(filter_results)
        filter_results("")

        layout.addWidget(table)
        layout.addWidget(label(
            self.tr_ui("読み込み元: ", "Source: ")
            + fasta_path, 10
        ))
        layout.addWidget(label(
            self.tr_ui(
                "変異位置はFASTA参照配列を基準とする1始まりの番号です。"
                "ProteinMPNNのスコアは実験的な活性や安定性を"
                "直接示すものではありません。",
                "Mutation positions are relative to the FASTA "
                "reference sequence (1-based). ProteinMPNN scores "
                "do not directly measure experimental activity "
                "or stability."
            ),
            11
        ))

        return widget

    def records_page(self):
        widget, layout = self.page(
            "研究記録",
            "再現性を確保するための実行環境と条件を確認します。"
        )
        if self.data:
            manifest = self.data["manifest"]
            layout.addWidget(label(
                "ProteinMPNN Commit: "
                + manifest.get("proteinmpnn_commit", "不明")
            ))
            layout.addWidget(label(
                "PyTorch: " + manifest.get("pytorch", "不明")
            ))
            layout.addWidget(label(
                "記録ファイル数: "
                + str(len(manifest.get("files", {})))
            ))
        else:
            layout.addWidget(label("研究記録を読み込めませんでした。"))
        layout.addStretch()
        return widget

    def settings_page(self):
        widget, layout = self.page(
            self.tr_ui("設定", "Settings"),
            self.tr_ui(
                "アプリと計算環境に関する情報です。",
                "Application and computing environment."
            )
        )
        layout.addWidget(label(self.tr_ui("操作画面：Windows / PySide6", "Interface: Windows / PySide6")))
        layout.addWidget(label(self.tr_ui("計算環境：WSL2 / Ubuntu", "Compute environment: WSL2 / Ubuntu")))
        layout.addWidget(label(self.tr_ui("現在のバージョン：0.1（閲覧専用）", "Version: 0.1 (read-only)")))
        layout.addStretch()
        return widget


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Studio()
    window.show()
    sys.exit(app.exec())


