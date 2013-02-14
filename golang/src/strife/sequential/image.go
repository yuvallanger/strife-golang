package sequential

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"io/ioutil"
)

func (model *Model) Save_snapshots_as_images() {
	for snapshot_i, snapshot := range model.Data_samples.Snapshots {
		img := board_strain_to_image(&snapshot.Data)
		imagefile, err := ioutil.TempFile(".", fmt.Sprintf("image-%v-%05d-png-", model.Start_Time, snapshot_i))
		if err != nil {
			panic(err)
		}
		png.Encode(imagefile, img)
		imagefile.Close()
	}
}

func board_strain_to_image(board_strain *Board_strain) (img *image.RGBA) {
	img = image.NewRGBA(image.Rect(0,
		0,
		len(*board_strain),
		len(*board_strain)))
	for row := range *board_strain {
		for col, strain := range (*board_strain)[row] {
			img.Set(row, col, color.RGBA{
				uint8(255 * r4strain[strain]),
				uint8(255 * s4strain[strain]),
				uint8(255 * g4strain[strain]),
				255})
		}
	}
	return img
}
