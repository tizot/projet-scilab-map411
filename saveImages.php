<?php
$matrix = array_diff(scandir('.'), array('.', '..'));

foreach ($matrix as $file) {
    if (!preg_match('#txt$#is', $file)) {
        continue;
    }

    $m = file_get_contents($file);
    $m = explode("\n", $m);
    $lines = count($m);
    $columns = count(explode(' ', $m[0]));

    $image = imagecreatetruecolor($columns, $lines);
    $gray = array();
    for ($i = 0; $i <= 255; $i++) {
        $gray[$i] = imagecolorallocate($image, $i, $i, $i);
    }

    $y = 0;
    foreach ($m as $l) {
        $y++;
        $x = 0;

        $line = explode(' ', $l);
        foreach ($line as $c) {
            $x++;

            $color = $gray[floor($c*255)];

            imagesetpixel($image, $x, $y, $color);
        }
        echo '.';
    }
    imagepng($image, $file . '.png');
    echo "\n";
}
