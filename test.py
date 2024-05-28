from PIL import Image

# 打开 PNG 图片
png_image = Image.open('icon.png')

# 转换并保存为 ICO
png_image.save('icon.ico')