package sequential

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"io/ioutil"
)

func (model *Model) SaveSnapshotsAsImages() {
	for snapshot_i, snapshot := range model.DataSamples.Snapshots {
		img := board_strain_to_image(&snapshot.Data)
		imagefile, err := ioutil.TempFile(".", fmt.Sprintf("image-%v-%05d-png-", model.StartTime, snapshot_i))
		if err != nil {
			panic(err)
		}
		png.Encode(imagefile, img)
		imagefile.Close()
	}
}

func board_strain_to_image(boardStrain *BoardStrain) (img *image.RGBA) {
	img = image.NewRGBA(image.Rect(0,
		0,
		len(*boardStrain),
		len((*boardStrain)[0])))
	for row := range *boardStrain {
		for col, strain := range (*boardStrain)[row] {
			img.Set(row, col, color.RGBA{
				uint8(255 * r4strain[strain]),
				uint8(255 * s4strain[strain]),
				uint8(255 * g4strain[strain]),
				255})
		}
	}
	return img
}
