from pathlib import Path

EXAMPLES_DIR = Path("examples")
OUTPUT_FILE = Path("examples.html")
EXTENSION = ".semisim"  # Change to your desired extension

BEGIN_TEMPLATE = """
<html>

<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Examples</title>
  <style>
    body {
      margin: 0;
      background: #f9f9f9;
    }

    .container {
      font-family: Arial, sans-serif;
      margin: 28px;
    }

    h1,
    h2 {
      color: #2c3e50;
    }

    .section {
      margin-bottom: 40px;
    }

    .example-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, 300px);
      gap: 20px;
    }

    .card {
      background: white;
      border: 1px solid #ddd;
      border-radius: 10px;
      box-shadow: 0 2px 5px rgba(0, 0, 0, 0.1);
      overflow: hidden;
      text-align: center;
      transition: transform 0.2s;
    }

    .card:hover {
      transform: scale(1.02);
    }

    .card img {
      width: 100%;
      height: 140px;
      object-fit: contain;
      background: #000;
    }

    .card a {
      display: block;
      padding: 10px;
      text-decoration: none;
      color: #2980b9;
      font-weight: bold;
    }
  </style>
  <link rel="stylesheet" href="style.css">
</head>

<body>
  <div class="container">
"""

FOLDER_TEMPLATE = """
    <div class="section">
      <h2>%s</h2>
      <div class="example-grid">
"""

THUMBNAIL_TEMPLATE = """        <div class="card">
          <img src="%s">"""

FILE_TEMPLATE = """          <a href="" onclick='app.load("%s");'>%s</a>
        </div>"""

FOLDER_END_TEMPLATE = """
      </div>
    </div>
"""

END_TEMPLATE = """
  </div>
</body>
</html>
"""

html = []

html.append(BEGIN_TEMPLATE)

preferred_order = [
    "Electrostatics",
    "Circuits",
    "Waves",
    "Semiconductor physics",
    "Diodes",
    "BJTs",
    "MOSFETs",
    "JFETs",
    "Misc semiconductors",
    "Digital logic",
    "Test",
]

subdirs = {d.name: d for d in EXAMPLES_DIR.iterdir() if d.is_dir()}
ordered = [subdirs.pop(name) for name in preferred_order if name in subdirs]
ordered.extend(sorted(subdirs.values(), key=lambda d: d.name))

for subdir in ordered:
    if not subdir.is_dir():
        continue

    html.append(FOLDER_TEMPLATE % subdir.name)

    for file in sorted(subdir.glob(f"*{EXTENSION}")):
        file_path = file.as_posix()
        thumbnail_path = file_path.replace("examples", "images/thumbnails", 1).replace(EXTENSION, ".png", 1)

        html.append(THUMBNAIL_TEMPLATE % thumbnail_path)
        html.append(FILE_TEMPLATE % (file_path, file.stem))

    html.append(FOLDER_END_TEMPLATE)

html.append(END_TEMPLATE)

with OUTPUT_FILE.open("w", encoding="utf-8") as f:
    f.write("\n".join(html))

print(f"Wrote {OUTPUT_FILE}")