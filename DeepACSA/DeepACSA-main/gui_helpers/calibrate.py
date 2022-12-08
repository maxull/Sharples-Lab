"""Python module which provides functions to calibrate US images."""

import math
import numpy as np
import cv2

mlocs = []

def region_of_interest(img, vertices):
    """Defines region of interest where ridges are searched.

    Arguments:
        Processed image containing edges,
        numpy array of regions of interest vertices / coordinates.

    Returns:
        Masked image of ROI containing only ridges detected by preprocessing.

    Example:
        >>>region_of_interest(preprocessed_image,
        np.array([(0,1), (0,2), (4,2), (4,7)], np.int32),)
    """
    mask = np.zeros_like(img)
    # channel_count = img.shape[2]
    match_mask_color = 255
    cv2.fillPoly(mask, vertices, match_mask_color)
    masked_image = cv2.bitwise_and(img, mask)
    return masked_image

def mclick(event, x_val, y_val, flags, param):
    """Detect mouse clicks for purpose of image calibration.

    Arguments:

    Returns:
        List of y coordinates of clicked points.
    """
    global mlocs

    # if the left mouse button was clicked, record the (x, y) coordinates
    if event == cv2.EVENT_LBUTTONDOWN:
        mlocs.append(y_val)
        mlocs.append(x_val)

def draw_the_lines(img, line):
    """Draws lines along the detected ridges

    Arguments:
        Original image,
        numpy array of detected lines by Hough-Transform.

    Returns:
        Original images containing lines.

    Example:
        >>>draw_the_lines(Image1.tif,
        np.array([[0 738 200 539]))
    """
    img = np.copy(img)
    # Creating empty image to draw lines on
    blank_image = np.zeros((img.shape[0], img.shape[1], 3), dtype=np.uint8)

    for x_1, y_1, x_2, y_2 in line:
        cv2.line(blank_image, (x_1, y_1), (x_2, y_2), (0, 255, 0),
                thickness=3)

    # Overlay image with lines on original images (only needed for plotting)
    img = cv2.addWeighted(img, 0.8, blank_image, 1, 0.0)
    return img

def calibrate_distance_efov(path_to_image: str, arg_muscle: str):
    """Calculates scalingline length of image based computed
        length of detected rigdes.

        Arguments:
            Path to image that should be analyzed.

        Returns:
            Length of scalingline (pixel).

        Example:
            >>>calibrate_distance_efov(C:/Desktop/Test/Img1.tif)
            571
    """
    image = cv2.imread(path_to_image)
    # Transform BGR Image to RGB
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    height = image.shape[0]
    width = image.shape[1]
    # Define ROI with scaling lines
    region_of_interest_vertices = [
        (10, height),
        (10, height*0.1),
        (width, height*0.1),
        (width, height)
    ]
    # Transform RGB to greyscale for edge detection
    gray_image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
    # Edge detecition
    canny_image = cv2.Canny(gray_image, 400, 600)
    cropped_image = region_of_interest(canny_image,
                                       np.array([region_of_interest_vertices],
                                                np.int32),)

    # For RF
    muscle = arg_muscle
    if muscle == "RF":
        lines = cv2.HoughLinesP(cropped_image,
                                rho=1,
                                theta=np.pi/180,
                                threshold=50,
                                lines=np.array([]),
                                minLineLength=325,
                                maxLineGap=3)
        if lines is None:
            return None, None
        # draw lines on image
        image_with_lines = draw_the_lines(image, lines[0])

    # For VL
    if muscle == "VL":
        lines = cv2.HoughLinesP(cropped_image,
                                rho=1,  # Distance of pixels in accumulator
                                theta=np.pi/180,  # Angle resolution
                                threshold=50,  # Only lines with higher vote
                                lines=np.array([]),
                                minLineLength=175,
                                maxLineGap=3)  # Gap between lines
        if lines is None:
            return None, None
        # draw lines on image
        image_with_lines = draw_the_lines(image, lines[0])

    # For GM / GL
    if muscle == "GL":
        lines = cv2.HoughLinesP(cropped_image,
                                rho=1,  # Distance of pixels in accumulator
                                theta=np.pi / 180,  # Angle resolution
                                threshold=50,  # Only lines with higher vote
                                lines=np.array([]),
                                minLineLength=250,
                                maxLineGap=5)
        if lines is None:
            return None, None
        # draw scaling lines on image
        image_with_lines = draw_the_lines(image, lines[0])

    else:
        lines = cv2.HoughLinesP(cropped_image,
                                rho=1,  # Distance of pixels in accumulator
                                theta=np.pi / 180,  # Angle resolution
                                threshold=50,  # Only lines with higher vote
                                lines=np.array([]),
                                minLineLength=250,
                                maxLineGap=5)
        if lines is None:
            return None, None
        # draw scaling lines on image
        image_with_lines = draw_the_lines(image, lines[0])

    # Calculate length of the scaling line
    scalingline = lines[0][0]
    point1 = [scalingline[0], scalingline[1]]
    point2 = [scalingline[2], scalingline[3]]
    scalingline_length = math.sqrt(((point1[0] - point2[0])**2)
                                   + ((point1[1] - point2[1])**2))

    return scalingline_length, image_with_lines

def calibrate_distance_static(nonflipped_img, spacing: str):
    """Calculates scalingline length of image based computed
        distance between two points on image and image depth.

    Arguments:
        Original(nonflipped) image with scaling lines on right border,
        distance between scaling points (mm).

    Returns:
        Length of scaling line (pixel).

    Example:
        >>>calibrate_distance_manually(Image, 5, 4.5, 0)
        5 mm corresponds to 95 pixels
    """
    # calibrate according to scale at the right border of image
    img2 = np.uint8(nonflipped_img)
    height = img2.shape[0]
    width = img2.shape[1]
    imgscale = img2[int(height*0.4):(height), (width-int(width*0.15)):width]

    # search for rows with white pixels, calculate median of distance
    calib_dist = np.max(np.diff(np.argwhere(imgscale.max(axis=1) > 150),
                                axis=0))

    if int(calib_dist) < 1:
        return None, None, None

    # calculate calib_dist for 10mm
    if spacing == "5":
        calib_dist = calib_dist * 2
    if spacing == "15":
        calib_dist = calib_dist * (2/3)
    if spacing == "20":
        calib_dist = calib_dist / 2

    #scalingline_length = depth * calib_dist
    scale_statement = '10 mm corresponds to ' + str(calib_dist) + ' pixels'

    return calib_dist, imgscale, scale_statement

def calibrate_distance_manually(nonflipped_img, spacing):
    """Calculates scalingline length of image based on manual specified
        distance between two points on image and image depth.

    Arguments:
        Original(nonflipped) image,
        distance between scaling points (mm).

    Returns:
        Length of scaling line (pixel).

    Example:
        >>>calibrate_distance_manually(Image, 5)
        5 mm corresponds to 99 pixels
    """
    img2 = np.uint8(nonflipped_img)

    # display the image and wait for a keypress
    cv2.imshow("image", img2)
    cv2.setMouseCallback("image", mclick)
    key = cv2.waitKey(0)

    # if the 'q' key is pressed, break from the loop
    if key == ord("q"):
        cv2.destroyAllWindows()

    global mlocs

    calib_dist = np.abs(math.sqrt((mlocs[3] - mlocs[1])**2 + (mlocs[2] - mlocs[0])**2))
    mlocs = []
    # calculate calib_dist for 10mm
    if spacing == 5:
        calib_dist = calib_dist * 2
    if spacing == 15:
        calib_dist = calib_dist * (2/3)
    if spacing == 20:
        calib_dist = calib_dist / 2

    # print(str(spacing) + ' mm corresponds to ' + str(calib_dist) + ' pixels')

    return calib_dist
