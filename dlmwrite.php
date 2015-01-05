<?php
if (!isset($argv[1])) {
    echo "Use: /Applications/MAMP/bin/php/php5.5.10/bin/php dlmwrite.php Basketball/frame10.png > img1.txt\n";
    exit();
}

$path = $argv[1];

$img = imagecreatefrompng($path);

$dimensions = getimagesize($path);
$width = $dimensions[0];
$height = $dimensions[1];

for ($y = 1; $y <= $height; $y++) {
    for ($x = 1; $x <= $width; $x++) {
        $rgb = imagecolorat($img, $x, $y);
        echo round($rgb/256, 3);
        if ($x < $width) {
            echo " ,";
        }
    }
    echo "\n";
}
