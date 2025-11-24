# =============================================================================
# Interactive Image Region Selector and Cropper
# 交互式图像区域选择和裁剪工具
# =============================================================================
# Purpose / 用途:
#   Manually select and crop square regions from images using mouse clicks
#   使用鼠标点击手动选择和裁剪图像的正方形区域
# 
# Usage Instructions / 使用说明:
#   1. Modify IMAGE_NAME below to select your target image / 修改下面的IMAGE_NAME选择目标图像
#   2. Run this script / 运行此脚本
#   3. Click two points to define a square region / 点击两个点定义正方形区域
#   4. The cropped region will be saved with consistent naming / 裁剪区域将以一致的命名保存
#   5. Press ESC to close the window and move to next image / 按ESC键关闭窗口并移至下一张图
#   6. Script will exit after processing the image / 脚本将在处理完图像后退出
# =============================================================================

import cv2
import numpy as np
import os

# =============================================================================
# Configuration / 配置
# =============================================================================
# Modify this variable to process different images / 修改此变量以处理不同的图像
# Example values: "400_1", "400_2", "400_3", "600_1", "600_2", "600_3", "800_1", "800_2", "800_3"
IMAGE_NAME = "400_1"

# Input directory / 输入目录
INPUT_DIR = "1_gray_colored_png"

# Output directory for cropped regions / 裁剪区域输出目录
OUTPUT_DIR = "2_catched"

# =============================================================================
# Auto-generated paths / 自动生成的路径
# =============================================================================
# Input image path / 输入图像路径
INPUT_IMAGE_PATH = os.path.join(INPUT_DIR, f"gray{IMAGE_NAME}.png")

# Output image path for cropped region / 裁剪区域输出路径
OUTPUT_ROI_PATH = os.path.join(OUTPUT_DIR, f"liti{IMAGE_NAME}.png")

# Output image path for processed region (with transparency) / 处理后区域输出路径（带透明度）
OUTPUT_NEW_ROI_PATH = os.path.join(OUTPUT_DIR, f"new_liti{IMAGE_NAME}.png")

print("=" * 80)
print("Configuration / 配置:")
print(f"  Image name / 图像名称: {IMAGE_NAME}")
print(f"  Input image / 输入图像: {INPUT_IMAGE_PATH}")
print(f"  Output ROI / 输出ROI: {OUTPUT_ROI_PATH}")
print(f"  Output processed / 输出处理后: {OUTPUT_NEW_ROI_PATH}")
print("=" * 80)

# =============================================================================
# Global Variables / 全局变量
# =============================================================================
# Store clicked coordinates / 存储点击的坐标
coords = []

# Create output directory if it doesn't exist / 如果输出目录不存在则创建
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")
    print(f"已创建输出目录: {OUTPUT_DIR}")

# =============================================================================
# Mouse Callback Function / 鼠标回调函数
# =============================================================================
def draw_line(event, x, y, flags, param):
    """
    Mouse callback function to handle click events and draw selection rectangle.
    鼠标回调函数，处理点击事件并绘制选择矩形。
    
    Parameters / 参数:
    ----------------
    event : int
        OpenCV mouse event type / OpenCV鼠标事件类型
    x, y : int
        Mouse cursor coordinates / 鼠标光标坐标
    flags : int
        OpenCV event flags / OpenCV事件标志
    param : any
        Additional parameters (not used) / 附加参数（未使用）
    
    Processing Logic / 处理逻辑:
    --------------------------
    1. On left mouse click, store the coordinate / 左键点击时存储坐标
    2. After two clicks, calculate square dimensions / 两次点击后计算正方形尺寸
    3. Draw rectangle on image / 在图像上绘制矩形
    4. Crop and save the region / 裁剪并保存区域
    5. Clear coordinates for next selection / 清空坐标以便下次选择
    
    WARNING / 警告:
    -------------
    DO NOT modify the cropping logic or coordinate calculations!
    不要修改裁剪逻辑或坐标计算！
    """
    global coords, img

    # Step 1: Capture left mouse button click / 步骤1：捕获左键点击
    if event == cv2.EVENT_LBUTTONDOWN:
        coords.append((x, y))

    # Step 2: Process after two clicks / 步骤2：两次点击后处理
    if len(coords) == 2:
        # Calculate square side length (use max of width/height) / 计算正方形边长（使用宽度/高度的最大值）
        # WARNING: DO NOT MODIFY! / 警告：不要修改！
        side = max(abs(coords[1][0] - coords[0][0]), abs(coords[1][1] - coords[0][1]))

        # Draw rectangle on image (green color, 1px thickness) / 在图像上绘制矩形（绿色，1像素粗细）
        cv2.rectangle(img, coords[0], (coords[0][0] + side, coords[0][1] + side), (0, 255, 0), 1)

        # Crop region of interest (ROI) / 裁剪感兴趣区域（ROI）
        # WARNING: DO NOT MODIFY the +1 and -1 offsets! / 警告：不要修改+1和-1偏移量！
        roi = img[coords[0][1]+1:coords[0][1] + side - 1, coords[0][0]+1:coords[0][0] + side - 1]

        # Save cropped region / 保存裁剪区域
        cv2.imwrite(OUTPUT_ROI_PATH, roi)
        print(f"\nROI saved / ROI已保存: {OUTPUT_ROI_PATH}")

        # Clear coordinates for next selection / 清空坐标以便下次选择
        coords.clear()

# =============================================================================
# Step 1: Interactive Selection Process / 步骤1：交互选择流程
# =============================================================================

# Load image / 加载图像
img = cv2.imread(INPUT_IMAGE_PATH)

if img is None:
    print(f"ERROR: Cannot load image / 错误：无法加载图像: {INPUT_IMAGE_PATH}")
    print("Please check if the file exists / 请检查文件是否存在")
    exit(1)
else:
    print(f"Image loaded successfully / 图像加载成功: {INPUT_IMAGE_PATH}")

# Create window / 创建窗口
cv2.namedWindow('image')

# Set mouse callback function / 设置鼠标回调函数
cv2.setMouseCallback('image', draw_line)

# Main display loop / 主显示循环
while True:
    # Display image / 显示图像
    cv2.imshow('image', img)
    
    # Wait for ESC key (27) to exit / 等待ESC键（27）退出
    if cv2.waitKey(20) & 0xFF == 27:
        break

# Close all windows / 关闭所有窗口
cv2.destroyAllWindows()

# =============================================================================
# Step 2: Apply Threshold-Based Transparency to Cropped Image
# 步骤2：对裁剪图像应用基于阈值的透明度
# =============================================================================

def catch(image_path, gray_level):
    """
    Apply threshold-based transparency to an image.
    对图像应用基于阈值的透明度处理。
    
    Parameters / 参数:
    ----------------
    image_path : str
        Path to the input image file / 输入图像文件的路径
    gray_level : int
        Grayscale threshold value (0-255) / 灰度阈值（0-255）
        Pixels below this value will become transparent / 低于此值的像素将变为透明
    
    Processing Logic / 处理逻辑:
    --------------------------
    1. Load image (supports color or grayscale) / 加载图像（支持彩色或灰度）
    2. Convert to grayscale if necessary / 如需要则转换为灰度图
    3. Create 4-channel image (RGBA) / 创建4通道图像（RGBA）
    4. Copy grayscale values to RGB channels / 将灰度值复制到RGB通道
    5. Set alpha channel to 0 for pixels below threshold / 将低于阈值的像素的alpha通道设为0
    6. Save result / 保存结果
    
    WARNING / 警告:
    -------------
    DO NOT modify the gray_level threshold or processing logic!
    不要修改gray_level阈值或处理逻辑！
    The threshold value of 64 is calibrated for specific visualization requirements.
    阈值64是为特定可视化需求校准的。
    
    Output / 输出:
    ------------
    4-channel PNG image with transparency applied / 应用了透明度的4通道PNG图像
    
    Example / 示例:
    -------------
    Input:  liti400_1.png (grayscale image)
    Output: new_liti400_1.png (RGBA image with transparency)
    
    Threshold Behavior / 阈值行为:
    ----------------------------
    - Pixels with value >= gray_level: Opaque (alpha = 255) / 不透明（alpha = 255）
    - Pixels with value < gray_level: Transparent (alpha = 0) / 透明（alpha = 0）
    """
    
    # Step 1: Load image / 步骤1：加载图像
    img = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)

    # Step 2: Convert to grayscale if needed / 步骤2：如需要则转换为灰度图
    if len(img.shape) == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Step 3: Create 4-channel image (RGBA) / 步骤3：创建4通道图像（RGBA）
    new_img = np.ones((img.shape[0], img.shape[1], 4), dtype=np.uint8) * 255

    # Step 4: Copy grayscale to RGB channels / 步骤4：将灰度值复制到RGB通道
    new_img[:, :, :3] = img[:, :, np.newaxis]

    # Step 5: Apply threshold-based transparency / 步骤5：应用基于阈值的透明度
    # WARNING: DO NOT MODIFY! / 警告：不要修改！
    new_img[img < gray_level, 3] = 0

    # Step 6: Save result / 步骤6：保存结果
    cv2.imwrite(OUTPUT_NEW_ROI_PATH, new_img)

# Check if ROI file exists / 检查ROI文件是否存在
if os.path.exists(OUTPUT_ROI_PATH):
    catch(OUTPUT_ROI_PATH, 64)
    print(f"Processed image saved / 处理后图像已保存: {OUTPUT_NEW_ROI_PATH}")
    print("\n" + "=" * 80)
    print("Processing complete! / 处理完成！")
    print("You can now run this script again with a different IMAGE_NAME")
    print("您现在可以使用不同的IMAGE_NAME再次运行此脚本")
    print("=" * 80)
else:
    print(f"ERROR: ROI file not found / 错误：ROI文件未找到: {OUTPUT_ROI_PATH}")
    print("Please complete the interactive selection first / 请先完成交互式选择")

