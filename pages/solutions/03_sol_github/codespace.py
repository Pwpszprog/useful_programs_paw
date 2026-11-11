import os
import urllib.request
from PIL import Image, ImageDraw, ImageFont

# 1. Pobieranie japońskiej czcionki (Noto Sans JP), jeśli nie istnieje w katalogu projektu
font_path = "NotoSansJP.ttf"
font_url = "https://github.com/google/fonts/raw/main/ofl/notosansjp/NotoSansJP%5Bwght%5D.ttf"

if not os.path.exists(font_path):
    print("Pobieranie czcionki japońskiej...")
    urllib.request.urlretrieve(font_url, font_path)
    print("Pobrano czcionkę NotoSansJP.ttf")

# 2. Tworzenie obrazka (800x400 px, białe tło)
width, height = 800, 400
image = Image.new("RGB", (width, height), color="white")
draw = ImageDraw.Draw(image)

# 3. Właściwy tekst: "dziękuję" (a-ri-ga-to-u)
text = "ありがとう"

# 4. Wczytanie pobranej czcionki
font = ImageFont.truetype(font_path, size=70)

# 5. Wyśrodkowanie tekstu
bbox = draw.textbbox((0, 0), text, font=font)
text_width = bbox[2] - bbox[0]
text_height = bbox[3] - bbox[1]
x = (width - text_width) // 2
y = (height - text_height) // 2 - bbox[1]

# 6. Rysowanie tekstu i zapis do pliku
draw.text((x, y), text, fill="#1e293b", font=font)
image.save("arigatou.png")
print("Zapisano plik jako arigatou.png")