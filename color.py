import dataclasses
import copy
import random
import matplotlib.colors as mcolors
import random
from PIL import ImageColor


@dataclasses.dataclass
class Color:
    color_name: str = ""
    hex: str = ""
    rgb: tuple = ()
    target: str = ""

"""
旧カラーコード
    "Europe": "#FF1234",
    "North_America": "#00FF00",
    "Oceania": "#0000FF",
    "Asia": "#00cccc",
    "Africa": "#ff00ff",
    "South_America": "#ffcc00",
"""

targets = {
    "continent": {
        "Europe": "red",
        "North_America": "lime",
        "Oceania": "blue",
        #         "Asia": "sorbus",
        "Asia": "aqua",
        "Africa": "fuchsia",
        #         "Africa": "violets are blue",
        "South_America": "yellow",
    },
    "continent4": {
        "Europe": "red",
        "North_America": "lime",
        "Oceania": "blue",
        "Asia": "aqua",
    },
    "Pango lineage": {
        "B.1.617.2": "fuchsia",
        "AY.4": "aqua",
        "AY.12": "blueviolet",
        "AY.3": "dark orange",
        "AY.6": "teal",
        "AY.9": "pink",
        "AY.5": "blue",
        "AY.7": "brown",
        "AY.11": "dark green",
    },  # AY3,1m AY.2,AY.10,AY.1,AY.8除外
    "month": {
        "11": "fuchsia",
        "12": "aqua",
        "01": "blueviolet",
        "02": "dark orange",
        "03": "teal",
        "04": "pink",
        "05": "blue",
        # "06": "",
        # "07": "brown",
        # "08": "dark green",
        # "09": "portage",
        # "10": "purple",
        # "11": "yellow",
    },
    "clade": {
        "G": "fuchsia",
        "GH": "aqua",
        "GR": "blueviolet",
        "GV": "dark orange",
        "GRY": "teal",
        "L": "pink",
        "O": "blue",
        "S": "brown",
        "V": "dark green",
    },
}


def create_color_code(n):
    targets = {}
    random.seed(0)
    while n > 0:
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        if f"0x{hex(r)[2:]}{hex(g)[2:]}{hex(b)[2:]}" not in targets:
            targets[f"0x{hex(r)[2:]}{hex(g)[2:]}{hex(b)[2:]}"] = [r, g, b]
            n -= 1
    targets["white"] = [255, 255, 255]
    targets["black"] = [0, 0, 0]
    return targets


def define_rgb_color(classes, target=None, seed=0):
    random.seed(seed)
    colors = copy.copy(mcolors.CSS4_COLORS)
    del colors["black"], colors["white"]

    colors = [Color(color_name=i, hex=v, rgb=ImageColor.getcolor(v, "RGB")) for i, v in colors.items()]

    if target is None:
        colors = random.sample(colors, k=len(classes))
        rgb_color = {}
        for i in range(len(classes)):
            colors[i].target = classes[i]
            rgb_color[classes[i]] = colors[i]
    else:
        rgb_color = {}
        for color in colors:
            for con, color_name in  targets[target].items():
                if color_name == color.color_name:
                    rgb_color[con] = color
    return rgb_color


# rgb_color = {
#     "torch red": [255, 18, 52],
#     "lime": [0, 255, 0],
#     "blue": [0, 0, 255],
#     "sorbus": [237, 107, 53],
#     "fuchsia": [255, 0, 255],
#     "yellow": [255, 255, 0],
#     "brass": [181, 166, 66],
#     "green": [0, 102, 0],
#     "kelly green": [51, 204, 0],
#     "aqua": [0, 255, 255],
#     "portage": [153, 153, 255],
#     "dark orange": [255, 127, 0],
#     "orange red": [255, 69, 0],
#     "sea green": [46, 139, 87],
#     "bakers chocolate": [92, 51, 23],
#     "barberry": [217, 217, 25],
#     "dark green": [0, 100, 0],
#     "indigo": [75, 0, 130],
#     "deeppink": [255, 20, 147],
#     "gray": [128, 128, 128],
#     "purple": [128, 0, 128],
#     "teal": [0, 128, 128],
#     "sandy brown": [244, 164, 96],
#     "brown": [143, 101, 82],
#     "darkviolet": [148, 0, 211],
#     "blueviolet": [138, 43, 226],
#     "pink": [255, 192, 203],
#     "violets are blue": [131, 102, 244],  # 9966FF
#     "robin egg blue": [0, 204, 204],  # 00cccc
#     "white": [255, 255, 255],
#     "black": [0, 0, 0],
# }

